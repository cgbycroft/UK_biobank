#!/bin/bash
# The name of the job, can be anything, simply used when displaying the list of running jobs
#$ -N UKBB-pca-UKBio-White_British
# Giving the name of the output log file
#$ -o Logs
#$ -j y
# One needs to tell the queue system to use the current directory as the working directory
# Or else the script may fail as it will execute in your top level home directory
#$ -cwd
#$ -V
#$ -P donnelly.prja -q short.qa


## A small program to compute projections given pre-computed loadings
flashproj=/well/ukbiobank/qcoutput/Src/coolibah/flashproj/flashproj

basedir=/well/ukbiobank/qcoutput.V2_QCed.sample-QC

## SNP loadings computed for release2 pipeline testing pipeline_v0.sh
## %%%%% check that this is correct
SNPloads=$basedir/data/PCA/b1__b11-b001__b095-autosome-sampleqc-fastpca-highquality-White_British_rare.snpload.map

## running in cluster
batchList=$basedir/QC-Scripts/batchList.txt

Batch=`head -n $SGE_TASK_ID $batchList | tail -n 1`

## or
## Choose coolibah:
#BatchNo=094
#i=1
#Batch=`head -n $i $batchList | tail -n 1`

Datadir=$basedir/data/ByBatch  # %%% <==== change these directories to point to batch-based SNP-QC'd plink base files
PCAdir=$basedir/data/PCA

## Usage: flashproj --bfile PlinkData --snploads LoadingsFile
## Options:
##   --freqA2  If speficied then FRQ = FRQ(A2). Otherwise FRQ = FRQ(A1).


## From the plink2 website: In the .frq file produced by --freq:
## A1    Allele 1 (usually minor)
## A2    Allele 2 (usually major)
## MAF   Allele 1 frequency
## Of course, since we always use the --keep-allele-freq option,
## MAF is the frequency of the reference allele


$flashproj --bfile $Datadir/$Batch-autosome-sampleqc --snploads $SNPloads --out $PCAdir/$Batch-PCA


## Append the sample IDs to create a final version of the PCA projections output file
paste <(cut -d' ' -f1 $Datadir/$Batch-autosome-sampleqc.fam) $PCAdir/$Batch-PCA.proj > $PCAdir/$Batch-PCA.pcs


## NOTE: Run plot-pca-UKBio.R in order to keep a unique version of the pcs combined across batches.
## The rest of the analysis to correct the heterozygosity for population structure is done in R.
