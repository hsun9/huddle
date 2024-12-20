"""
  Hua Sun
  12/20/24 v0.4
  
"""

import glob
import sys
import pandas as pd
import os
import argparse



parser = argparse.ArgumentParser()
parser.add_argument('-d', '--dir', type=str, default='.', help='cellranger outs')
parser.add_argument('-o', '--outdir', type=str, default='out_summary', help='out dir')

args = parser.parse_args()



def Main():
    func = args.func
    if not os.path.isdir(args.outdir):
            os.mkdir(args.outdir)

    # summary per sample
    samples = get_first_level_folders(args.dir)
    for x in samples:
        print(x)
        f_csv = f'{args.dir}/{x}/outs/summary.csv'
        QC_MultiomeSummaryFile(f_csv, f'{args.outdir}/{x}.summary.log')

    # merge
    MergeSummaryFile(args.outdir)



"""
  Set Func.
"""

def get_first_level_folders(path):
    folders = []
    for entry in os.scandir(path):
        if entry.is_dir():
            folders.append(entry.name)
    return folders



## Cellragner-arc count 'summary.csv'
def QC_MultiomeSummaryFile(f_data, outfile):
  df = pd.read_csv(f_data, header=None)
  df = df.T
  df.columns = df.iloc[0]
  df.drop(df.index[0], inplace=True)
  df.set_index('Sample ID', inplace=True)

  name = df.columns[0]
  report = pd.DataFrame(
      [['Estimated number of cells','500-10000','','FAIL'],

      ['Feature linkages detected','-','','-'],
      ['Linked genes','-','','-'],
      ['Linked peaks','-','','-'],

      ['ATAC Confidently mapped read pairs','>0.8','','FAIL'],
      ['ATAC Fraction of genome in peaks','<0.75','','FAIL'],
      ['ATAC Fraction of high-quality fragments in cells','>0.4','','FAIL'],

      ['ATAC Fraction of high-quality fragments overlapping peaks','>0.25','','FAIL'],
      ['ATAC Fraction of transposition events in peaks in cells','>0.25','','FAIL'],
      ['ATAC Mean raw read pairs per cell','>5000','','FAIL'],
      ['ATAC Median high-quality fragments per cell','>100','','FAIL'],
      ['ATAC Non-nuclear read pairs','<0.1','','FAIL'],
      
      ['ATAC TSS enrichment score','>5','','FAIL'],
      ['ATAC Valid barcodes','>0.85','','FAIL'],
      
      ['ATAC Number of peaks','-','','-'],
      ['ATAC Percent duplicates','-','','-'],
      ['ATAC Q30 bases in barcode','-','','-'],
      ['ATAC Q30 bases in read 1','-','','-'],
      ['ATAC Q30 bases in read 2','-','','-'],
      ['ATAC Q30 bases in sample index i1','-','','-'],
      ['ATAC Sequenced read pairs','-','','-'],
      ['ATAC Unmapped read pairs','-','','-'],
      ['ATAC Fraction of high-quality fragments overlapping TSS','-','','-'],

      ['GEX Fraction of transcriptomic reads in cells','>0.6','','FAIL'],
      ['GEX Mean raw reads per cell','>5000','','FAIL'],
      ['GEX Median UMI counts per cell','>100','','FAIL'],
      ['GEX Median genes per cell','>500','','FAIL'],

      ['GEX Reads mapped antisense to gene','<0.3','','FAIL'],
      ['GEX Reads mapped confidently to genome','>0.5','','FAIL'],
      ['GEX Reads mapped confidently to intergenic regions','<0.3','','FAIL'],
      ['GEX Reads mapped confidently to transcriptome','>0.5','','FAIL'],
      ['GEX Reads mapped to genome','>0.8','','FAIL'],
      ['GEX Reads with TSO','<0.25','','FAIL'],

      ['GEX Valid barcodes','>0.8','','FAIL'],

      ['GEX Percent duplicates','-','','-'],
      ['GEX Q30 bases in UMI','-','','-'],
      ['GEX Q30 bases in barcode','-','','-'],
      ['GEX Q30 bases in read 2','-','','-'],
      ['GEX Q30 bases in sample index i1','-','','-'],
      ['GEX Q30 bases in sample index i2','-','','-'],
      ['GEX Reads mapped confidently to exonic regions','-','','-'],
      ['GEX Reads mapped confidently to intronic regions','-','','-'],
      ['GEX Sequenced read pairs','-','','-'],
      ['GEX Total genes detected','-','','-'],
      ['GEX Valid UMIs','-','','-']],

      columns=['Catalog', '10xStandard', name, f'{name}.qc'])

  for catalog in report.Catalog:
    val = df.loc[catalog, name]
    report.loc[report.Catalog==catalog,name] = df.loc[catalog, name]

    val = float(val)

    if catalog == 'Estimated number of cells':
      if (val > 500) & (val < 10000):
        report.loc[report.Catalog==catalog,f'{name}.qc'] = 'PASS'
    if catalog == 'ATAC Confidently mapped read pairs':
      if val > 0.8:
        report.loc[report.Catalog==catalog,f'{name}.qc'] = 'PASS'
    if catalog == 'ATAC Fraction of genome in peaks':
      if val < 0.75:
        report.loc[report.Catalog==catalog,f'{name}.qc'] = 'PASS'
    if catalog == 'ATAC Fraction of high-quality fragments in cells':
      if val > 0.4:
        report.loc[report.Catalog==catalog,f'{name}.qc'] = 'PASS'
    if catalog == 'ATAC Fraction of high-quality fragments overlapping peaks':
      if val > 0.25:
        report.loc[report.Catalog==catalog,f'{name}.qc'] = 'PASS'
    if catalog == 'ATAC Fraction of transposition events in peaks in cells':
      if val > 0.25:
        report.loc[report.Catalog==catalog,f'{name}.qc'] = 'PASS'
    if catalog == 'ATAC Mean raw read pairs per cell':
      if val > 5000:
        report.loc[report.Catalog==catalog,f'{name}.qc'] = 'PASS'
    if catalog == 'ATAC Median high-quality fragments per cell':
      if val > 100:
        report.loc[report.Catalog==catalog,f'{name}.qc'] = 'PASS'
    if catalog == 'ATAC Non-nuclear read pairs':
      if val < 0.1:
        report.loc[report.Catalog==catalog,f'{name}.qc'] = 'PASS'
    if catalog == 'ATAC TSS enrichment score':
      if val > 5:
        report.loc[report.Catalog==catalog,f'{name}.qc'] = 'PASS'
    if catalog == 'ATAC Valid barcodes':
      if val > 0.85:
        report.loc[report.Catalog==catalog,f'{name}.qc'] = 'PASS'
    if catalog == 'GEX Fraction of transcriptomic reads in cells':
      if val > 0.6:
        report.loc[report.Catalog==catalog,f'{name}.qc'] = 'PASS'
    if catalog == 'GEX Mean raw reads per cell':
      if val > 5000:
        report.loc[report.Catalog==catalog,f'{name}.qc'] = 'PASS'
    if catalog == 'GEX Median UMI counts per cell':
      if val > 100:
        report.loc[report.Catalog==catalog,f'{name}.qc'] = 'PASS'
    if catalog == 'GEX Median genes per cell':
      if val > 500:
        report.loc[report.Catalog==catalog,f'{name}.qc'] = 'PASS'
    if catalog == 'GEX Reads mapped antisense to gene':
      if val < 0.3:
        report.loc[report.Catalog==catalog,f'{name}.qc'] = 'PASS'
    if catalog == 'GEX Reads mapped confidently to intergenic regions':
      if val < 0.3:
        report.loc[report.Catalog==catalog,f'{name}.qc'] = 'PASS'
    if catalog == 'GEX Reads mapped confidently to transcriptome':
      if val > 0.5:
        report.loc[report.Catalog==catalog,f'{name}.qc'] = 'PASS'
    if catalog == 'GEX Reads mapped to genome':
      if val > 0.8:
        report.loc[report.Catalog==catalog,f'{name}.qc'] = 'PASS'
    if catalog == 'GEX Reads mapped confidently to genome':
      if val > 0.5:
        report.loc[report.Catalog==catalog,f'{name}.qc'] = 'PASS'
    if catalog == 'GEX Reads with TSO':
      if val < 0.25:
        report.loc[report.Catalog==catalog,f'{name}.qc'] = 'PASS'
    if catalog == 'GEX Valid barcodes':
      if val > 0.8:
        report.loc[report.Catalog==catalog,f'{name}.qc'] = 'PASS'

  #print(report)
  #n_pass = report.loc[report.iloc[:,3]=='PASS'].shape[0]/report.shape[0] * 100
  n_pass = report.loc[report.iloc[:,3]=='PASS'].shape[0]/report.loc[report.iloc[:,3].isin(['PASS','FAIL'])].shape[0] * 100
  n_pass = round(n_pass, 2)
  print(f'Percentage of PASS: {n_pass}%')
  report = report._append(dict(zip(report.columns,['Percentage of PASS', '.', '.', n_pass])), ignore_index=True)
  report.to_csv(outfile, sep='\t', index=False)





# merge summary files
def MergeSummaryFile(dir_qc):
    files = glob.glob(f'{dir_qc}/*/*.summary.log')
    i = 0
    df_qc = []
    for f in files:
        d_temp = []
        if i == 0:
            d_temp = pd.read_csv(f, sep='\t')
            df_qc = d_temp
            i = 1
        else:
            d_temp = pd.read_csv(f, sep='\t', usecols=[0,2,3])
            df_qc = df_qc.merge(d_temp, on='Catalog')

        print(f'[INFO] file:{f} - {d_temp.shape}')

    df_qc.to_csv(f'{dir_qc}/merged_all.qc_summary.xls', sep='\t', index=False)






if __name__ == '__main__':
  Main()






