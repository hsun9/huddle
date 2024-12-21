# Hua Sun
# 2024-12-20 v0.3

RUN='local'

## getOptions
while getopts "C:R:G:N:O:L:F:n:" opt; do
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
    L)
      LIBRARIES=$OPTARG
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

if [ ! -f ${LIBRARIES} ]; then
    echo "[ERROR] The input librar ${LIBRARIES} does not exist !"
    exit 1
fi

# set reference
reference=''
if [[ $GENOME == 'hg38' ]]; then
    reference=${REF_10x_GRCh38}
elif [[ $GENOME == 'mm10' ]]; then
    reference=${REF_10x_MM10}
else
    echo '[ERROR] No matched species ! Please set -S human or mouse '
    exit 1
fi


cd ${OUTDIR}

# ATAC+GeneExp
# local mode: stable
if [[ $RUN == 'local' ]]; then
  ${CELLRANGER_ARC} count --id=${NAME} \
                          --reference=${reference} \
                          --libraries=${LIBRARIES} \
                          --localcores=32 \
                          --localmem=64
fi



if [[ $RUN == 'local-fc' ]]; then
    ${CELLRANGER_ARC} count --id=${NAME} \
                      --reference=${reference} \
                      --libraries=${LIBRARIES} \
                      --force-cells=${ncell} \
                      --localcores=16 \
                      --localmem=64
fi




# cluster mode: unstable. it depends on cluster server plots
if [[ $RUN == 'lsf' ]]; then
  ${CELLRANGER_ARC} count --id=${NAME} \
                          --reference=${reference} \
                          --libraries=${LIBRARIES} \
                          --jobmode=lsf \
                          --mempercore=8
fi



