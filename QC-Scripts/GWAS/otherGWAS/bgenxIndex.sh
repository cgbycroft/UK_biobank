#!/bin/bash
#$ -cwd
#$ -N bgenxIndex_hrc+uk10k.I4.v1.1
#$ -o Logs
#$ -j y
#$ -V

#########
# Create indexed bgen files in v1.1 format
#########
bgenix=/apps/well/bgenix/20160708/bgenix
basedir=/well/ukbiobank/qcoutput.V2_QCed.sample-QC

#mkdir $basedir/data/imputed
#mkdir $basedir/data/imputed/Logs

cd $basedir/data/imputed

chr=$SGE_TASK_ID

echo "Creating indexed bgen for chromosome $chr"

inputBgen=/well/ukbiobank/imputation/final/full/bgen/hrc+uk10k_sorted_v1.1/chr${chr}.hrc+uk10k.I4.v1.1

# make symbolic link in my directory
ln -s $inputBgen.bgen $basedir/data/imputed
ln -s $inputBgen.sample $basedir/data/imputed


echo $bgenix $basedir/data/imputed/chr${chr}.hrc+uk10k.I4.v1.1.bgen
time $bgenix $basedir/data/imputed/chr${chr}.hrc+uk10k.I4.v1.1.bgen



# TO RUN FROM $basedir/data/imputed : qsub -P donnelly.prjb -q long.qb -t 1-22 bgenxIndex.sh
