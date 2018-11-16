#!/bin/bash
#$ -cwd
#$ -N bolt_lreg_v1
#$ -o Logs
#$ -j y
#$ -pe shmem 2
#$ -V
#$ -t 1-22

# Downloaded from https://data.broadinstitute.org/alkesgroup/BOLT-LMM/downloads/ on May 12th 2016
bolt=/well/donnelly/ukbiobank_project_8874/clare/src/BOLT-LMM_v2.2/bolt 
basedir=/well/ukbiobank/qcoutput.V2_QCed.sample-QC

vers=v1 # Oxford QC
pheno=Standing.height  # <==== This should exist in the phenoFile

chrom=$SGE_TASK_ID

mkdir $basedir/data/GWAS/otherGWAS/$pheno
mkdir $basedir/data/GWAS/otherGWAS/$pheno/BOLTLMM.$vers

genoData=$basedir/data/Combined/byChrom/b1__b11-b001__b095-autosome-oxfordqc-chr$chrom
phenoFile=$basedir/data/GWAS/otherGWAS/PhenotypesForBOLT-$vers.txt

# phenoFile and snpSet was created using subsetPhenotypes.R
snpSet=$basedir/data/PCA/b1__b11-b001__b095-autosome-sampleqc-fastpca-init-highquality.prune.in  # just use SNPs from original PCA :) NOT ACTUALLY USED IN LREG!!

outFile=$basedir/data/GWAS/otherGWAS/$pheno/BOLTLMM.$vers/$pheno-BOLT-LREG-$vers-$chrom.out

# WHOLE GENOME MEMORY NEEDED: 74,100,452,428 bytes = 5 C-NODES  <==== Sample QC SNPs
# WHOLE GENOME MEMORY NEEDED: 96,739,226,940  = 7 C-NODES   <==== Oxford QC SNPs

$bolt \
--bfile=$genoData \
--phenoFile=$phenoFile \
--phenoCol=$pheno \
--covarFile=$phenoFile \
--covarCol=Array \
--covarCol=Inferred.Gender \
--qCovarCol=PC{1:8} \
--numThreads=2 \
--statsFile=$outFile \
--covarMaxLevels=111 \
--geneticMapFile=/well/donnelly/ukbiobank_project_8874/clare/src/BOLT-LMM_v2.2/tables/genetic_map_hg19.txt.gz \
--verboseStats

# qsub -q short.qc -P donnelly.prjc runBOLT_v1-LREG.sh

