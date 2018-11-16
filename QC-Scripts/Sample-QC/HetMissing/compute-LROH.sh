#!/bin/bash
# The name of the job, can be anything, simply used when displaying the list of running jobs
#$ -N UKB-hetmissing
# Giving the name of the output log file
#$ -o Logs
#$ -j y
#$ -cwd
#$ -V
#$ -P donnelly.prjc -q short.qc


plink=/well/ukbiobank/qcoutput/Src/plink-v1.90b3k-x86_64/plink

########### The variable basedir is set in the pipeline_vx.sh file
inputdir=$basedir/data/Combined
outputdir=$basedir/data/Relatedness
###########

minDistance=1000
echo $inputdir

time $plink --bfile $inputdir/b1__b11-b001__b095-autosome-sampleqc --homozyg-kb $minDistance --out $outputdir/b1__b11-b001__b095-autosome-sampleqc-$minDistance.KB
