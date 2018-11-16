#!/bin/bash
#$ -cwd
#$ -N UKB-ibd-plink
#$ -o Logs
#$ -j y
#$ -V
#$ -P donnelly.prjc
#$ -q short.qc
#$ -t 1-50

basedir=/well/ukbiobank/qcoutput.V2_QCed.sample-QC
plink=/well/ukbiobank/qcoutput/Src/plink-v1.90b3k-x86_64/plink

######
threshold=0.08838835 # equivalent to KING's 3rd degree threshold (i.e 2*0.04419417 where kinship coefficient = 1/2^(9/2))
######


# Use the same set of SNPs as with KING!
input=$basedir/data/Combined/b1__b11-b001__b095-autosome-sampleqc-filtered.0.003

output=$basedir/data/Relatedness/b1__b11-b001__b095-autosome-sampleqc-filtered.0.003-White_British

# just compute IBD on white british subset, and only those that ended up in the pruned kinship table! These are selected using select-IBD-samples.R
samples=$basedir/QC-Scripts/Sample-QC/Relatedness/plink-IBD-samples.txt



## $plink --bfile $input --keep $samples --genome full --min $threshold --parallel $SGE_TASK_ID --out $output

echo $plink --bfile $input --keep $samples --genome full --parallel $SGE_TASK_ID 50 --min $threshold --out $output

time $plink --bfile $input --keep $samples --genome full --parallel $SGE_TASK_ID 50 --min $threshold --out $output
