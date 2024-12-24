# Hua Sun
# v2.4  12/23/24
# based on v2.3.2


library(Seurat)
library(Signac)
library(dplyr)
library(stringr)
library(stringi)
library(harmony)
library(GetoptLong)
library(this.path)
path <- dirname(this.path())
source(paste0(path,'/utils.R'))

set.seed(47)



rdir <- 'out_cellranger_arc'
ref <- 'hg38'
min_cells <- 3
outdir <- 'out_multiome_integrated_h5'

GetoptLong(
    "rdir=s",       "dir",
    "ref=s",       "ref",
    "min_cells=i", "min_cells",
    "outdir=s",    "outdir"
)


dir.create(outdir)


# set
macs2 <- '/research/groups/mackgrp/home/common/Software/miniconda3/envs/macs2/bin/macs2'
annotation <- Annotations(ref)


multiome_list <- list()
folder_names <- list.dirs(rdir, full.names=F, recursive=F)

for (name in folder_names){
    print(name)
    h5 <- paste0(rdir, '/', name, '/outs/filtered_feature_bc_matrix.h5')
    fragment <- paste0(rdir, '/', name, '/outs/atac_fragments.tsv.gz')

    seu <- CreateObjectFrome10xH5(h5=h5, min_cells=min_cells, atac_fragment=fragment, ref=ref, name=name, annotation=annotation)

    seu$Sample <- name
    seu$ref <- ref

    # make data list
    multiome_list <- append(multiome_list, seu)
}


print(multiome_list)
saveRDS(multiome_list, file=paste0(outdir,'/rawdata.list.rds'))


multiome_filtered_list <- lapply(X = multiome_list, FUN = function(x) { x <- FilterCells_Multiome(x) })
saveRDS(multiome_filtered_list, file=paste0(outdir,'/filtered.list.rds'))


multiome_filtered_list <- lapply(X = multiome_list, FUN = function(x) { x <- FilterDoublet(x) })
saveRDS(multiome_filtered_list, file=paste0(outdir,'/filtered.doublet.list.rds'))



##---------------- snRNA: 2-1.Merge ----------------##
# https://github.com/satijalab/seurat/issues/4896
# https://satijalab.org/seurat/articles/sctransform_v2_vignette
# return.only.var.genes set default 
# perform SCTransform normalization
snrna_list <- base::lapply(X = multiome_filtered_list, FUN = function(x) { x <- SCTransform(x, vst.flavor = "v2", variable.features.n = 3000) })

# sct-integration & pca
features <- SelectIntegrationFeatures(object.list = snrna_list, nfeatures = 3000)
write.table(features, file = paste0(outdir, '/snrna_selectedIntegrationFeatures.xls'), sep = "\t", quote=F, col.names = NA)


# Merge normalized samples for harmony
samples <- c()
for (x in snrna_list){ samples <- append(samples, unique(x$Sample)) }

# https://satijalab.org/seurat/articles/merge_vignette
combined_snrna <- merge(x=snrna_list[[1]], y=snrna_list[2:length(snrna_list)], merge.data = TRUE, add.cell.ids = samples)

# save merged
saveRDS(combined_snrna, file=paste0(outdir,'/snrna_merged.rds'))




##---------------- snRNA: 2-2.Harmony ----------------##

DefaultAssay(combined_snrna) <- "SCT"
VariableFeatures(combined_snrna) <- features

harmony_rna <- combined_snrna %>%
                RunPCA(assay = "SCT", npcs = 50) %>%
                RunHarmony(group.by.vars="Sample", reduction = "pca", assay.use = "SCT", reduction.save="harmony_rna")

# save harmony data before umap
saveRDS(harmony_rna, file=paste0(outdir,'/snrna_harmony.pre.rds'))


# default seed.use=42
# https://satijalab.org/seurat/articles/sctransform_v2_vignette.html
harmony_rna <- RunUMAP(harmony_rna, reduction = "harmony_rna", 
                        reduction.name = "umap.rna", reduction.key = "rnaUMAP_", assay = "SCT", dims = 1:30) %>%
                FindNeighbors(reduction = "harmony_rna", dims = 1:30) %>%
                FindClusters(resolution=0.3)

# save harmony data
saveRDS(harmony_rna, file=paste0(outdir,'/snrna_harmony.rds'))
write.table(harmony_rna@meta.data, file = paste0(outdir, '/snrna_harmony.metadata.xls'), sep = "\t", quote=F, col.names = NA)

# finial barcodes
filtered_code <- stri_replace_last_fixed(row.names(harmony_rna@meta.data), '_', '#')
writeLines(filtered_code, paste0(outdir, '/filtered.barcode.archr_format.out'))


# plot umap
p <- DimPlot(harmony_rna, group.by = 'Sample', reduction = "umap.rna", pt.size = 0.1)
pdf(paste0(outdir, '/snrna_harmony.umap.pdf'), width = 5.5, height = 4, useDingbats=FALSE)
print(p)
dev.off()



# remove processed
snrna_list <- NULL
combined_snrna <- NULL





##---------------- snATAC: 3-1.Call peaks using MACS2 ----------------##
snatac_macs2_list <- lapply(X = multiome_filtered_list, FUN = function(x) { x <- CallPeaksUsingMACS2(x, macs2, annotation) })
saveRDS(snatac_macs2_list, file=paste0(outdir, '/snatac_macs2.list.rds'))





##---------------- snATAC: 3-2.Combine peaks & recall fragment ----------------##
# from macs2 results 
atac_peak <- c()
for (obj in snatac_macs2_list){
    name <- unique(obj$Sample)

    # make peak bed for multiome data
    peaks <- as.data.frame(obj@assays$peaks@ranges)[,1:3]
    colnames(peaks) <- c('chr', 'start', 'end')
    
    # convert to genomic ranges
    gr <- makeGRangesFromDataFrame(as.data.frame(peaks))

    atac_peak[[name]] <- gr
}

print(names(atac_peak))



# 2.Create a unified set of peaks to quantify in each dataset
combined.peaks <- GenomicRanges::reduce(x = unlist(GRangesList(atac_peak)))

# 3.Filter out bad peaks based on length
peakwidths <- width(combined.peaks)
combined.peaks <- combined.peaks[peakwidths < 10000 & peakwidths > 20]

snatac_peak_list <- lapply(X = snatac_macs2_list, FUN = function(x) { x <- MakeNewATACObject(x, combined.peaks, annotation) })

saveRDS(snatac_peak_list, file=paste0(outdir, '/snatac_peaks.list.rds'))




##---------------- snATAC: 3-2b. Add gene activity ----------------##
pfm <- CallPFM_JASPAR2022()
saveRDS(pfm, file = paste0(outdir, "/pfm.jaspar2022.", ref, ".rds"))

print('snatac_peakExtra_list')
snatac_peakExtra_list <- lapply(X = snatac_peak_list, FUN = function(x) { x <- AddGeneAct_Motif(x, ref, pfm) })
saveRDS(snatac_peakExtra_list, file=paste0(outdir, '/snatac_peaksExtra.list.rds'))




##---------------- snATAC: 3-3.merge snATAC ----------------##

#samples <- c()
#for (x in snatac_peakExtra_list){ samples <- append(samples, unique(x$Sample)) }
combined_snatac <- merge(x=snatac_peakExtra_list[[1]], y=snatac_peakExtra_list[2:length(snatac_peakExtra_list)], add.cell.ids = samples)
combined_snatac[["peaks"]]
saveRDS(combined_snatac, file=paste0(outdir, '/snatac_merged.rds'))

# remove processed
snatac_macs2_list <- NULL
snatac_peak_list <- NULL
snatac_peakExtra_list <- NULL




##---------------- snATAC: 3-4. harmony snATAC ----------------##
# https://stuartlab.org/signac/articles/merging.html
# dims = 2:50

DefaultAssay(combined_snatac)<-"peaks"

combined_snatac <- combined_snatac %>% 
            RunTFIDF() %>%
            FindTopFeatures(min.cutoff = 50) %>%
            RunSVD() %>%
            RunUMAP(dims = 2:50, reduction = "lsi", reduction.name = "umap.atac_merge", reduction.key="atacMergedUMAP_")


p <- DimPlot(combined_snatac, group.by = 'Sample', reduction = "umap.atac_merge", pt.size = 0.1)
pdf(paste0(outdir, '/snatac_merged.umap.pdf'), width = 5.5, height = 4, useDingbats=FALSE)
print(p)
dev.off()



# 6.batch correction using harmony
harmony_atac <- combined_snatac %>%
                RunHarmony(group.by.vars = 'Sample', reduction = 'lsi', assay.use = 'peaks', , reduction.save="harmony_atac", project.dim = FALSE ) %>%
                RunUMAP(reduction = 'harmony_atac', reduction.name = "umap.atac", reduction.key = "atacUMAP_", dims = 2:30) %>%
                FindNeighbors(reduction = 'harmony_atac', dims = 2:30) %>%
                FindClusters(resolution = 0.3)


# save harmony data
saveRDS(harmony_atac, file=paste0(outdir,'/snatac_harmony.rds'))
write.table(harmony_atac@meta.data, file = paste0(outdir, '/snatac_harmony.metadata.xls'), sep = "\t", quote=F, col.names = NA)

p <- DimPlot(harmony_atac, group.by = 'Sample', reduction = "umap.atac", pt.size = 0.1)
pdf(paste0(outdir, '/snatac_harmony.umap.pdf'), width = 6.5, height = 5, useDingbats=FALSE)
print(p)
dev.off()



# remove processed
combined_snatac <- NULL




##---------------- Combind snRNA and snATAC ----------------##
# https://stuartlab.org/signac/articles/pbmc_multiomic.html
# https://github.com/satijalab/seurat/issues/6094

# Combine RNA and ATAC
integrated_multiome <- harmony_atac
integrated_multiome[['RNA']] <- harmony_rna[['RNA']]
integrated_multiome[['SCT']] <- harmony_rna[['SCT']]
integrated_multiome[['pca']] <- harmony_rna[['pca']]
integrated_multiome[['harmony_rna']] <- harmony_rna[['harmony_rna']]
integrated_multiome[['umap.rna']] <- harmony_rna[['umap.rna']]

# Compute a joint neighbor graph
integrated_multiome <- FindMultiModalNeighbors(
                object = integrated_multiome, 
                reduction.list = list("harmony_rna", "harmony_atac"), 
                dims.list = list(1:30, 2:30), 
                modality.weight.name = "SCT.weight",
                verbose = TRUE
            )

integrated_multiome <- RunUMAP(
                object = integrated_multiome, 
                nn.name = "weighted.nn", 
                reduction.name = "wnn.umap", 
                reduction.key = "wnnUMAP_",
                assay = "SCT",
                seed.use = 47,
                verbose = TRUE
            )
    
# cluster
integrated_multiome <- FindClusters(integrated_multiome, graph.name="wsnn", resolution=0.4)

# save harmony data
saveRDS(integrated_multiome, file=paste0(outdir,'/multiome_integrated.rds'))
write.table(integrated_multiome@meta.data, file = paste0(outdir, '/multiome_integrated.metadata.xls'), sep = "\t", quote=F, col.names = NA)

p <- DimPlot(integrated_multiome, reduction = 'wnn.umap', group.by = 'Sample', pt.size = 0.1)
pdf(paste0(outdir, '/multiome_integrated.wnnUMAP.pdf'), width = 6.5, height = 5, useDingbats=FALSE)
print(p)
dev.off()







