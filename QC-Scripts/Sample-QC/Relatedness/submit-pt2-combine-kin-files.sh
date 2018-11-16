#!/bin/bash
# The name of the job, can be anything, simply used when displaying the list of running jobs
#$ -N UKB-kinship-combine-init
# Giving the name of the output log file
#$ -o Logs
#$ -j y
# One needs to tell the queue system to use the current directory as the working directory
# Or else the script may fail as it will execute in your top level home directory
#$ -cwd
#$ -V
#$ -P donnelly.prjc -q short.qc


## The input file is a list of the {nBatches choose 2} batches in the dataset:
## The first and the second columns specify the full path to the two batches,
## in plink binary format (all individuals, only SNPs for sample QC)
## The third column specifies the _directory_ to store the King results
batchPairs=b1__b11-b001__b095-pair_batches
PairsOfBatches=$batchPairs.txt

## Kinship coefficients are computed in pairs of batches on the cluster.
##                   This script puts the related pairs in one big file.


## Usage: merge-kinship-coefficients(ListPairsBatches,KinshipCoeffsOut,KinshipFilePrefix)
##     or merge-kinship-coefficients(ListPairsBatches,KinshipCoeffsOut,SampleOutliersIn,KinshipFilePrefix) where
##     the optional third argument is a list of sample IDs to exclude from any related pairs.
##     For example, this might be a list of all samples with self-declared mixed ancestry or
##     hetmiss (heterozygosity/missingness) outliers.
## CLARE: Will simply exclude these pairs later using filter-king.R

KinshipCoeffs=$basedir/data/Relatedness/$batchPairs-init.kin0
perl postking/merge-kinship-coefficients.pl $PairsOfBatches $KinshipCoeffs degree3-init.kin0
