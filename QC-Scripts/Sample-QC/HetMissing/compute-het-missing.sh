#!/bin/bash
# The name of the job, can be anything, simply used when displaying the list of running jobs
#$ -N UKB-hetmissing-filtered-snps
# Giving the name of the output log file
#$ -o Logs
#$ -j y
#$ -cwd
#$ -V
#$ -P donnelly.prjc -q short.qc


plink=/well/ukbiobank/qcoutput/Src/plink-v1.90b3k-x86_64/plink
basedir=/well/ukbiobank/qcoutput.V2_QCed.sample-QC

########### The variable basedir is set in the pipeline_vx.sh file
inputdir=$basedir/data/Combined
threshold=0.003 #this is the SNPload threshold for filtering SNPs in KING
###########

echo $inputdir

# Compute Het and missing just on snps used in KING (from PCA round 1)
SNPsToInclude=$basedir/QC-Scripts/PCA/pca-UKBio/b1__b11-b001__b095-autosome-sampleqc-fastpca-highquality-init-snpsToKeepPCs-$threshold.txt
# 93,512 SNPs

time $plink --bfile $inputdir/b1__b11-b001__b095-autosome-sampleqc --extract $SNPsToInclude --het gz --missing gz --make-bed --out $inputdir/b1__b11-b001__b095-autosome-sampleqc-filtered.$threshold

