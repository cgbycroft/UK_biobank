#!/bin/bash
#$ -cwd
#$ -N interim-compare-frequencies
#$ -o Logs
#$ -j y
#$ -V
#$ -P donnelly.prja
#$ -q short.qa
#$ -t 1-26

basedir=/well/ukbiobank/qcoutput.V2_QCed.sample-QC
plink=/well/ukbiobank/qcoutput/Src/plink-v1.90b3k-x86_64/plink

outDir=/well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/Combined

chroms=`cat <(seq 1 22) <(printf "X\nXY\nY\nMT\n")`
#chr=`cat <(seq 1 22) <(printf "X\nXY\nY\nMT\n") | head -n $SGE_TASK_ID | tail -n 1` 

for chr in $chroms;
do \
    
    GenotypesForReleaseDir=/well/ukbiobank/expt/V2_QCed.export/data/calls_export/v3
    GenotypesForReleaseFile=V2_QCed.export.UKBiLEVEAX_b1-b11.Batch_b001_b095.chr$chr


    echo $plink --bfile $GenotypesForReleaseDir/$GenotypesForReleaseFile --keep-allele-order --freqx --out $outDir/$GenotypesForReleaseFile

    time $plink --bfile $GenotypesForReleaseDir/$GenotypesForReleaseFile --keep-allele-order --freqx --out $outDir/$GenotypesForReleaseFile


done;
    
