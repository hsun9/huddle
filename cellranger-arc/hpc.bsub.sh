
MEM=$1;shift
CORES=$1;shift
NAME=$1;shift

DIR=`pwd`
LOG=$DIR/logs
mkdir -p $LOG

bsub -q large_mem -P ${NAME} -J ${NAME} -n ${CORES} -R "rusage[mem=${MEM}000] span[hosts=1]" -oo ${LOG}/${NAME}.log -eo ${LOG}/${NAME}.err "$@"
