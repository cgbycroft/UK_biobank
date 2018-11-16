#!/bin/bash
#$ -N file-indexes
# Giving the name of the output log file
#$ -o Logs
#$ -j y
#$ -cwd
#$ -V


basedir=/well/ukbiobank/qcoutput.V2_QCed.sample-QC
batchlist=/well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/batchList.txt

Batch=`head -n ${SGE_TASK_ID} $batchlist | tail -n 1`
log2File=/well/ukbiobank/source/V2_All/$Batch/AxiomGT1.cnv.log2ratio.final.txt.gz
bafFile=/well/ukbiobank/source/V2_All/$Batch/AxiomGT1.cnv.baf.final.txt.gz

outputdir=$basedir/data/CNV/$Batch

echo $Batch

mkdir $outputdir

python3 get-file-indexes.py --gzfile $log2File --out $outputdir/AxiomGT1.cnv.log2ratio.final

python3 get-file-indexes.py --gzfile $bafFile --out $outputdir/AxiomGT1.cnv.baf.final
