#!/bin/bash
#$ -cwd
#$ -N dupe-concordance
#$ -o Logs
#$ -j y
#$ -V

#######
# This script computes raw differences in genotypes between duplicates in the UK Biobank.
#######

basedir=/well/ukbiobank/qcoutput.V2_QCed.sample-QC
plink=/well/ukbiobank/qcoutput/Src/plink-v1.90b3k-x86_64/plink

################
dup=$SGE_TASK_ID
################

# The genotype data to use
genoData=$basedir/data/Combined/b1__b11-b001__b095-autosome-sampleqc

# Which version of the king output are we using?
kingRunPrefix=b1__b11-b001__b095-pair_batches.filtered-samples-pruned-200

# All the twins and duplicates (file generated by find-families.R)
duplicates=$basedir/QC-Scripts/Sample-QC/Relatedness/${kingRunPrefix}-duplicates-twins.txt



# Just compute distances for all of them and prune out genuine twins later.
cat <(tail -n +2 $duplicates | head -n $dup | tail -n 1 | cut -f2 -d' ') <(tail -n +2 $duplicates | head -n $dup | tail -n 1 | cut -f1 -d' ') | awk '{print $0,$0}' > $basedir/QC-Scripts/Sample-QC/Relatedness/keep.$dup.txt


# NOTE: Only compute on markers where both individuals are non-missing. Takes a couple of minutes to run for two individuals!

$plink --bfile $genoData --keep $basedir/QC-Scripts/Sample-QC/Relatedness/keep.$dup.txt --geno 0 --distance --out $basedir/data/Relatedness/${kingRunPrefix}-duplicates-genetic-distances-$dup

rm $basedir/data/Relatedness/${kingRunPrefix}-duplicates-genetic-distances-$dup.dist.id


# Write out R-readible file
## $plink --bfile $genoData --keep keep.$dup.txt --recode A-transpose --out test-$dup

# diff=584
# 3392 variants removed due to missing genotype data (--geno).
# 602484 variants and 2 people pass filters and QC.