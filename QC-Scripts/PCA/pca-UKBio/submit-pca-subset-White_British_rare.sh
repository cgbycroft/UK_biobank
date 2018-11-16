#!/bin/bash

plink=/well/ukbiobank/qcoutput/Src/plink-v1.90b3k-x86_64/plink
flashpca=/well/ukbiobank/qcoutput/Src/flashpca/flashpca
nThreads=4
basedir=/well/ukbiobank/qcoutput.V2_QCed.sample-QC

#######################
## I have picked one batch in order to test the script
## In general, use plink to merge many (all?) batches,
## at the autosomal SNPs only.
##
## I will also assume that all SNP QC has already been applied,
## It is also important to remove SNPs that are entirely missing
## in any batches (although with 100+ batches it might be okay
## id the SNP is missing in one batch? If you do decide to use
## such SNPs, it might happen that some PCs separate some batches
## from the rest.)
##
## I will also assume that relevant sample QC has already been applied,
## to exclude samples with high missingness, unusual heterozygosite, and
## closely related individuals (or at least, one of each pair of closely
## related individuals.)
Input=$basedir/data/Combined/b1__b11-b001__b095-autosome-sampleqc
Output=$basedir/data/PCA/b1__b11-b001__b095-autosome-sampleqc-fastpca-White_British_rare


## EXCLUDE SNPS
## A list of 23 SNP regions to exclude before principal component analysis
highLDregions=RegionsHighLD/highload.regions.txt

#maf=0.025  # originally this was 2.5%, which seems quite strict... try 1%
#maf=0.01  <===== This was the MAF used in the GWAS analyses etc.

####################
maf=0.001  # <===== This is just an experiment
####################

miss=0.015  # minimum missing rate for SNPs after sample filtering



## EXCLUDE SAMPLES
# (list created with pca-sample-filters.R). Related individuals etc.
Outliers=pca-UKBio-sample-exclusions-White_British.txt


## SUBSET PLINK DATASET
$plink --bfile $Input --remove $Outliers --maf $maf --geno $miss --freq --make-bed --keep-allele-order --out $Output-highquality


## PRUNE SNPs IN HIGH LD
Input=$Output-highquality


## I have chosen the parameters so that there are about 100,000 left after pruning
## I just repeated this part a few times until there are about that many SNPs in the prune.in list
$plink --bfile $Input --exclude range $highLDregions --indep-pairwise 1000 80 0.1 --out $Input
# v0:
# 0.1 --> ~200,000 SNPs
# 0.05 --> ~130,000 SNPs
# v1:
# 0.1 --> 111,455 SNPs  # with 77,448 unrelated samples
# 0.1 --> 106,198 SNPs
# As above but MAF set to 1% --> 146,991 SNPs


## Create a subset of the genotype data
$plink --bfile $Input --extract $Input.prune.in --make-bed --keep-allele-order --out $Input-pruned


#### ---> now run fastpca
