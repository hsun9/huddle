"""

  Hua Sun
  2024-12-20 v0.3

"""


import argparse
import os
import re
import pandas as pd

# collect input arguments
parser = argparse.ArgumentParser()
parser.add_argument('-v', '--ver', type=str, default='v9', help='cellranger version v7.2/v9')
parser.add_argument('-m', '--model', type=str, default='local', help='local/local-fc/lsf')
parser.add_argument('-nc', type=str, default=8000, help='--force-cells number')
parser.add_argument('-t', '--table', type=str, default='input.txt', required=True, help='data table for run multiple samples')

parser.add_argument('-o', '--outdir', type=str, default='.', help='out dir')

args = parser.parse_args()




def Main():
    path = os.path.dirname(__file__)
    f_config = f'{path}/config.{args.ver}.ini'
    src = f'{path}/src'
    f_bsub = f'{path}/src/hpc.bsub.sh'

    CheckDir('src', src)
    CheckFile('config', f_config)

    outdir = args.outdir
    if (outdir == '.'):
        outdir = f'out_cellranger.{args.ver}'

    if not os.path.isdir(outdir):
        os.mkdir(outdir)
    

    # set script
    script = f'{src}/cellranger_{args.ver}.count.sh'
    CheckFile('script', script)

    CallSeqCount_MultipleSamples(f_bsub, f_config, script, args.model, args.nc, args.table, outdir)    

    



"""
    Set Func.
"""

## Check file (absolute path)
def CheckFile(tag, f_path):
    if not os.path.isfile(f_path):
        print(f'[ERROR] The {tag} {f_path} does not exist !')
        exit()

## Check dir (absolute path)
def CheckDir(tag, dir_path):
    if not os.path.isdir(dir_path):
        print(f'[ERROR] The {tag} {dir_path} does not exist !')
        exit()


## Set script
def SetScript(pipeline, src, version):
    

    return script




## CallSeqCount_MultipleSamples
def CallSeqCount_MultipleSamples(f_bsub, f_config, script, model, ncell, f_table, outdir):

    info = pd.read_csv(f_table, sep='\t', names=['ref', 'name', 'fq_dir'])

    for index, row in info.iterrows():
        ref = row['ref']
        name = row['name']
        name = name.replace('.', '_')
        fq_dir = row['fq_dir']
        
        print(ref, name, fq_dir)

        cmd = f'sh {f_bsub} 64 16 scrna_{name} bash {script} -C {f_config} -R {model} -n {ncell} -G {ref} -F {fq_dir} -N {name} -O {outdir}'
        os.system(cmd)






if __name__ == '__main__':
    Main()



