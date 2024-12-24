# 11/1/23

# Set Library
library(Seurat)
library(Signac)
library(GenomicRanges)
library(future)

library(GenomeInfoDb)
library(EnsDb.Mmusculus.v79)
library(EnsDb.Hsapiens.v86)
library(dplyr)
library(ggplot2)
library(limma)
library(biovizBase)

library(BSgenome.Mmusculus.UCSC.mm10)
library(BSgenome.Hsapiens.UCSC.hg38)
library(SeuratDisk)
library(DoubletFinder)

library(chromVAR)
library(JASPAR2024)
library(TFBSTools)
library(motifmatchr)

library(stringr)
library(stringi)

library(patchwork)
library(data.table)

library(GetoptLong)

options(scipen = 999)


# set
macs2 <- '/research/groups/mackgrp/home/common/Software/miniconda3/envs/macs2/bin/macs2'

# set annotation genome
Annotations <- function(ref){
    ensdb <- EnsDb.Mmusculus.v79
    if (ref == 'hg38'){ ensdb <- EnsDb.Hsapiens.v86 }
    annotation <- GetGRangesFromEnsDb(ensdb = ensdb)
    seqlevels(annotation) <- paste0('chr', seqlevels(annotation))
    genome(annotation) <- ref

    annotation
}




# create seurat object from h5
CreateObjectFrome10xH5 <- function(h5=NULL, min_cells=min_cells, atac_fragment='atac_fragment.tsv.gz', ref='mm10', name='Multiome', annotation=NULL)
{
    # extract ATAC data
    print('1.read h5')
    seu <- Read10X_h5(h5)
    atac_counts <- seu$Peaks

    # extract RNA data
    rna_counts <- seu$`Gene Expression`
    # read rna seurat obj.
    #seu <- readRDS(rds)

    # mitochondria
    print('2.set mt')
    mt_pattern <- "^mt-"
    rb_pattern <- "^Rp[sl]"

    if (ref == 'hg38'){ 
        mt_pattern <- "^MT-" 
        rb_pattern <- "^RP[SL]"
    }
    print(mt_pattern)
    print(rb_pattern)

    print('3.creat rna obj')
    # a. RNA - Create Seurat object
    seu <- CreateSeuratObject(counts=rna_counts, assay="RNA", project=name, min.cells=min_cells, min.features=200)
    seu <- PercentageFeatureSet(seu, pattern = mt_pattern, col.name = "percent.mt")
    seu <- PercentageFeatureSet(seu, pattern = rb_pattern, col.name = "percent.ribo")

    print('4.read atac data')
    # b. Add in the ATAC-seq data
    # get gene annotations for mm10 and use peaks in standard chromosomes (chr1-chrM)
    # format to grange form to run macs2
    grange.counts <- StringToGRanges(rownames(atac_counts), sep = c(":", "-"))
    grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
    atac_counts <- atac_counts[as.vector(grange.use), ]


    print('5.creat atac obj')
    # add in the ATAC-seq data
    chromatin <- CreateChromatinAssay(
        counts = atac_counts,
        sep = c(":", "-"),
        genome = ref,
        fragments = atac_fragment,
        annotation = Annotations(ref)
    )

    # cell will use (match with snRNA)
    cell_rna <- rownames(seu@meta.data)
    cell_atac <- colnames(chromatin@counts)
    cell.use <- intersect(cell_rna, cell_atac)

    seu <- subset(seu, cells=cell.use)
    seu[["ATAC"]] <- subset(chromatin, cells=cell.use)
    print(paste('[INFO] Used cells -', length(cell.use)))

    ## snATAC
    DefaultAssay(seu) <- "ATAC"
    seu <- NucleosomeSignal(seu)
    seu <- TSSEnrichment(seu)
    #seu <- TSSEnrichment(seu, fast=FALSE) # slow, but allow plotting the accessibility profile at the TSS
    print('[INFO] Completed TSSEnrichment')

    # blacklist ratio
    if (ref == 'mm10'){
        seu$blacklist_ratio <- FractionCountsInRegion(
            object = seu, assay = 'ATAC',
            regions = blacklist_mm10
        )
    } 
       
    if (ref == 'hg38'){
        seu$blacklist_ratio <- FractionCountsInRegion(
            object = seu, assay = 'ATAC',
            regions = blacklist_hg38_unified
        )        
    }
    

    return(seu)
}




# raw data QC
QCPlotMultiome <- function(obj){
    p <- VlnPlot(
      object = obj,
      features = c('nCount_ATAC', 'nFeature_ATAC', 'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal', 'nCount_RNA', 'nFeature_RNA', 'percent.mt', 'percent.ribo'),
      pt.size = 0.01,
      ncol = 5
    ) & 
    theme(plot.title = element_text(size = 8, face = "bold"), axis.text=element_text(size = 8)) & 
    theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())

    p
}



# old
UpperVal_nCount <- function(values)
{
    Q <- quantile(values, probs = c(.1, .9), na.rm = T)
    Q1 <- Q[1]
    Q3 <- Q[2]
    IQR <- Q3 - Q1
    upper <- as.integer(Q3 + 1.5*IQR)

    return(upper)
}


UpperVal_nFeature <- function(values)
{
    Q <- quantile(values, probs = c(.2, .8), na.rm = T)
    Q1 <- Q[1]
    Q3 <- Q[2]
    IQR <- Q3 - Q1
    upper <- as.integer(Q3 + 1.5*IQR)

    return(upper)
}





# Filter cells
FilterCells_Multiome <- function(obj)
{   
    ref <- unique(obj$ref)
    print(ref)

    # min TSS
    #min_tss <- round(quantile(obj$TSS.enrichment, probs = .1, na.rm=T)[[1]])
    #if (min_tss < 2){ min_tss = 2}
    #print(min_tss)

    # max mt
    #max_mt <- round(quantile(obj$percent.mt, probs = .96, na.rm=T)[[1]])
    #if (max_mt > 10){ max_mt = 10}
    #if (max_mt == 0){ max_mt = 1}
    #print(max_mt)

    # max black list
    #max_bl <- round(quantile(obj$blacklist_ratio, probs = .96, na.rm=T)[[1]], 2)
    #if (max_bl > 0.05){ max_bl = 0.05}
    #print(max_bl)

    # max nCount_RNA ~50000
    #nc <- round(quantile(obj$nCount_RNA, probs = .996, na.rm=T)[[1]])
    #x <- 10^(nchar(nc)-1)
    #max_ncount_rna <- round(nc/x) * x
    #if (max_ncount_rna > 1){ max_ncount_rna = 100000 }
    #print(max_ncount_rna)

    obj[["log10GenesPerUMI_RNA"]] <- round(log10(obj$nFeature_RNA) / log10(obj$nCount_RNA), 2)
    
    # filter out low quality cells by QC plot
    obj <- subset(obj,
            subset = nCount_ATAC > 3000 &
            nCount_ATAC < 100000 &
            nucleosome_signal < 4 &
            TSS.enrichment > 3 &
            TSS.enrichment < 15 &
            blacklist_ratio < 0.03 &
            nCount_RNA > 2000 &
            nCount_RNA < 100000 &
            nFeature_RNA > 800 &
            nFeature_RNA < 8000 &
            percent.mt < 5 &
            log10GenesPerUMI_RNA > 0.8
        )

    print(dim(obj))
    obj
}



# Meta data log
MetadataLog <- function(obj){
    meta <- obj@meta.data

    catalog <- c('nCount_ATAC.min', 'nCount_ATAC.max', 
        'TSS.enrichment.min', 'TSS.enrichment.max', 'nucleosome_signal.max', 'blacklist_ratio.max',
        'nCount_RNA.min', 'nCount_RNA.max', 'nFeature_RNA.min', 'nFeature_RNA.max',
        'percent.mt.max', 'log10GenesPerUMI_RNA.min'
        )

    value <- c(min(meta$nCount_ATAC), max(meta$nCount_ATAC),
        min(meta$TSS.enrichment), max(meta$TSS.enrichment), max(meta$nucleosome_signal), max(meta$blacklist_ratio),
        min(meta$nCount_RNA), max(meta$nCount_RNA), min(meta$nFeature_RNA), max(meta$nFeature_RNA), 
        max(meta$percent.mt), min(meta$log10GenesPerUMI_RNA))

    d_log <- data.frame('catalog' = catalog, 'value' = as.numeric(value))

    #options(scipen = 999)
    d_log$value <- round(d_log$value, 2)

    d_log
}



# DoubletFinder
# https://rpubs.com/kenneditodd/doublet_finder_example
# https://raw.githubusercontent.com/ucdavis-bioinformatics-training/2022-July-Single-Cell-RNA-Seq-Analysis/main/data_analysis/scRNA_Workshop-PART7.Rmd
FilterDoublet <- function(obj)
{   
    DefaultAssay(obj) <- 'RNA'
    
    obj <- obj %>%
            NormalizeData() %>%
            FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
            ScaleData() %>%
            RunPCA()

    # Determine percent of variation associated with each PC
    pct <- obj[["pca"]]@stdev / sum(obj[["pca"]]@stdev) * 100
    # Calculate cumulative percents for each PC
    cumu <- cumsum(pct)
    # Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
    co1 <- which(cumu > 90 & pct < 5)[1]
    # Determine the difference between variation of PC and subsequent PC
    co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
    # Minimum of the two calculation
    min.pc <- min(co1, co2)

    obj <- obj %>%
            FindNeighbors(reduction="pca", dims = 1:min.pc) %>%
            FindClusters(resolution = 0.1) %>%
            RunUMAP(dims=1:min.pc)
          
    # 10x Single Cell 3' Gene Expression v3.1 assay
    # https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
    n_total_cell = nrow(obj@meta.data)
    doublet_rate = n_total_cell/1000*0.0075
    print(paste0('[INFO] Total cells: ', n_total_cell))
    print(paste0('[INFO] Estimated doublet rate: ', doublet_rate))

    # v2.0.4 'paramSweep' not 'paramSweep_v3'
    sweep.res <- DoubletFinder::paramSweep(obj, PCs = 1:min.pc, sct = FALSE)
    sweep.stats <- DoubletFinder::summarizeSweep(sweep.res, GT = FALSE)
    bcmvn <- DoubletFinder::find.pK(sweep.stats)
    pK <- bcmvn %>% dplyr::filter(BCmetric == max(BCmetric)) %>% dplyr::select(pK)
    pK <- as.numeric(as.character(pK[[1]]))

    # Homotypic Doublet Proportion Estimate 
    homotypic.prop <- DoubletFinder::modelHomotypic(obj@meta.data$seurat_clusters)
    nExp_poi <- round(doublet_rate*nrow(obj@meta.data))
    nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

    # Run DoubletFinder with varying classification stringencies 
    # v2.0.4 'doubletFinder' not 'doubletFinder_v3'
    obj <- DoubletFinder::doubletFinder(obj, PCs = 1:min.pc, pN = 0.25, pK = pK, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = FALSE)
    colnames(obj@meta.data) <- sub("^DF.classifications.*", "doubletfinder", colnames(obj@meta.data))

    # 
    sample_name <- unique(obj$orig.ident)
    pdf(paste0('doublet.', sample_name, '.umap.rna.pdf'), width = 6, height = 5, useDingbats=FALSE)
    p <- DimPlot(obj, group.by = 'doubletfinder')
    print(p)
    dev.off()

    # output for metadata with doublet results
    write.table(obj@meta.data, file = paste0('doublet.', sample_name, '.metadata.xls'), sep = "\t", quote=F, col.names = NA)

    # filter doublet from 
    cell.use <- colnames(subset(x = obj, subset = doubletfinder == 'Singlet'))
    obj <- subset(obj, cells=cell.use)

    return(obj)
}





# MACS2
CallPeaksUsingMACS2 <- function(obj=NULL, macs2='MACS2', annotation=NULL)
{
    DefaultAssay(obj) <- "ATAC"
    ref <- unique(obj$ref)

    effective_genome_size <- 2.3e+09

    if (ref == 'hg38'){ effective_genome_size <- 2.7e+09 }
    
    peaks <- CallPeaks(obj, macs2.path = macs2, effective.genome.size = effective_genome_size)

    # remove peaks on nonstandard chromosomes and in genomic blacklist regions
    peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")
    if (ref == 'mm10'){ peaks <- subsetByOverlaps(x = peaks, ranges = blacklist_mm10, invert = TRUE) }
    if (ref == 'hg38'){ peaks <- subsetByOverlaps(x = peaks, ranges = blacklist_hg38_unified, invert = TRUE) }

    # quantify counts in each peak
    macs2_counts <- FeatureMatrix(
        fragments = Fragments(obj),
        features = peaks,
        cells = colnames(obj)
    )

    obj[["peaks"]] <- CreateChromatinAssay(
        counts = macs2_counts,
        annotation = annotation
    )
   
    return(obj)
}


# JASPAR2022
CallPFM_JASPAR2022 <- function(){
    pfm <- getMatrixSet(x = JASPAR2022, opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE))
    pfm
}

# JASPAR2024
CallPFM_JASPAR2024 <- function(){
    pfm <- getMatrixSet(x = JASPAR2024, opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE))
    pfm
}


# re-call fragment
MakeNewATACObject <- function(obj=NULL, combined_peaks=NULL, annotation=NULL)
{
    # create new fragment obj
    sample <- unique(obj$Sample)
    metadata <- obj@meta.data
    barcode <- row.names(metadata)

    # make count based on combined peaks
    counts <- FeatureMatrix(fragments = Fragments(obj), features = combined_peaks, cells = barcode)
    
    # make object
    assay <- CreateChromatinAssay(counts, fragments=Fragments(obj), annotation = annotation)
    obj <- CreateSeuratObject(assay, assay = "peaks", meta.data=metadata)

    return(obj)
}




# add gene activity & motif information
AddGeneAct_Motif <- function(obj=NULL, ref='', pfm='')
{   
    print('GeneActivity')
    # https://stuartlab.org/signac/articles/pbmc_vignette.html
    # obj <- RunTFIDF(obj) %>% FindTopFeatures( min.cutoff="q5") %>% RunSVD()
    obj <- RunTFIDF(obj) %>% FindTopFeatures( min.cutoff="q0") %>% RunSVD()
    gene.activity <- GeneActivity(obj)

    print('CreateAssayObject')
    obj[['GA']] <- CreateAssayObject(counts = gene.activity)
    obj <- NormalizeData(object = obj, assay = 'GA', normalization.method = 'LogNormalize', scale.factor = median(obj$nCount_GA))


    # 2. motif
    # set genome
    print('[INFO] Call motif ...')
    genome <- BSgenome.Mmusculus.UCSC.mm10
    if (ref == 'hg38'){ genome <- BSgenome.Hsapiens.UCSC.hg38 }

    obj <- AddMotifs(object = obj, genome = genome, pfm = pfm)
    obj <- RunChromVAR(object = obj, genome = genome)

    return(obj)
}







