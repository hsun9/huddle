#!/bin/sh

# Hua Sun
# v0.2

## getOptions
RUN='local'

while getopts "C:R:G:N:O:F:n:" opt; do
  case $opt in
    C)
      CONFIG=$OPTARG
      ;;
    R)
      RUN=$OPTARG
      ;;
    G)
      GENOME=$OPTARG
      ;;
    N)
      NAME=$OPTARG
      ;;
    O)
      OUTDIR=$OPTARG
      ;;
    F)
      FQ_DIR=$OPTARG
      ;;
    n)
      ncell=$OPTARG
      ;;
    :)
      echo "Option -$OPTARG requires an argument." >&2
      exit 1
      ;;
  esac
done

source ${CONFIG}

mkdir -p $OUTDIR

if [ -d ${OUTDIR}/${NAME} ]; then
    echo "[ERROR] The name directory ${OUTDIR}/${NAME} exists !"
    exit 1
fi


# set reference
ref_gex=''
if [[ $GENOME == 'hg38' ]]; then
    ref_gex=${REF_GEX_10x_GRCh38}
elif [[ $GENOME == 'mm10' ]]; then
    ref_gex=${REF_GEX_10x_MM10}
else
    echo '[ERROR] No matched species ! Please set -S human or mouse '
    exit 1
fi


cd ${OUTDIR}


# local mode: stable
if [[ $RUN == 'local' ]]; then
  ${CELLRANGER} count --id=${NAME} \
                      --transcriptome=${ref_gex} \
                      --fastqs=${FQ_DIR} \
                      --localcores=16 \
                      --localmem=64
fi


if [[ $RUN == 'local-fc' ]]; then
    ${CELLRANGER} count --id=${NAME} \
                      --transcriptome=${ref_gex} \
                      --fastqs=${FQ_DIR} \
                      --force-cells=${ncell} \
                      --localcores=16 \
                      --localmem=64
fi


# cluster mode - but it depends on cluster
if [[ $RUN == 'lsf' ]]; then
  ${CELLRANGER} count --id=${NAME} \
                      --transcriptome=${ref_gex} \
                      --fastqs=${FQ_DIR} \
                      --jobmode=lsf \
                      --mempercore=8
fi

# GeneExp



