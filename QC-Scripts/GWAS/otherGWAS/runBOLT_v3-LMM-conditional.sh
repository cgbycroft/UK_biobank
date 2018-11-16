#!/bin/bash
#$ -cwd
#$ -N bolt_lmm_conditional_versionNumber
#$ -o Logs
#$ -j y
#$ -pe shmem 3
#$ -V

# Downloaded from https://data.broadinstitute.org/alkesgroup/BOLT-LMM/downloads/ on May 12th 2016
bolt=/well/donnelly/ukbiobank_project_8874/clare/src/BOLT-LMM_v2.2/bolt 
basedir=/well/ukbiobank/qcoutput.V2_QCed.sample-QC

vers=versionNumber
pheno=phenotypeForThisFile

### this sub-run is:
# BOLT version: LMM default - let BOLT decide on effect sizes.

mkdir $basedir/data/GWAS/otherGWAS/$pheno
mkdir $basedir/data/GWAS/otherGWAS/$pheno/BOLTLMM.$vers

genoDataAuto=$basedir/data/Combined/b1__b11-b001__b095-autosome-oxfordqc
genoDataSex=$basedir/data/Combined/b1__b11-b001__b095-sexchrom-oxfordqc

phenoFile=$basedir/data/GWAS/otherGWAS/PhenotypesForBOLT-$vers.txt  # using White british subset. Note: some individuals will be excluded in the genotype data

# phenoFile and snpSet was created using subsetPhenotypes.R
snpSet=$basedir/data/PCA/b1__b11-b001__b095-autosome-sampleqc-fastpca-White_British-highquality.prune.in  # just use SNPs from original PCA :) NOT ACTUALLY USED IN LREG!!  Note also that this does not include any SNPs on the sex chromosomes.

######## set up conditional SNPs etc. 
# regionsFile is created in region-plots.R and then conditional-analysis-setup.R
regionsFile=$pheno-BOLT-LMM-$vers.out.chrgenome.maf0.001.miss0.05.pruned-recombRegions-snpRegionList.txt

#SGE_TASK_ID=66

chrom=`cut -f1 -d' ' $regionsFile | tail -n +2 | head -n $SGE_TASK_ID | tail -n 1`
region=`cut -f2 -d' ' $regionsFile | tail -n +2 | head -n $SGE_TASK_ID | tail -n 1`
start=`cut -f3 -d' ' $regionsFile | tail -n +2 | head -n $SGE_TASK_ID | tail -n 1`
end=`cut -f4 -d' ' $regionsFile | tail -n +2 | head -n $SGE_TASK_ID | tail -n 1`

# set of SNPs to exclude: keep only SNPs in the region specified in $regionsFile
snpConditionFile=`cut -f5 -d' ' $regionsFile | tail -n +2 | head -n $SGE_TASK_ID | tail -n 1`
nsnps=`wc -l $snpConditionFile | cut -f1 -d' '`

awk -v s="$start" -v e="$end" -v c=$chrom ' ( $4 < s ) || ( $4>e ) || ( $1!=c ) {print $2}' $genoDataAuto.bim | grep -wvF -f $snpSet  > snp.exclusion.$chrom.$region.txt


# create temporary phenotype file
genotypesConditionalFile=genotypes-topVariants-$vers-phenoOrder-chr$chrom.txt

head -n 1 $genotypesConditionalFile | tr ' ' '\n' | grep -w -n -f $snpConditionFile > hitRegions/conditional.snp.order.$chrom.$region.txt

cols=`cut -f1 -d: hitRegions/conditional.snp.order.$chrom.$region.txt`
touch header$chrom.$region
o=1
for c in $cols; do echo SNP$o >> header$chrom.$region; o=`expr $o + 1`; done;
cols2=`echo $cols | tr ' ' ','`
header=`cat header$chrom.$region | tr '\n' ' '`
rm header$chrom.$region

paste -d' ' \
<(cat $phenoFile ) \
<(cat <(echo $header) <(cut -f$cols2 -d' ' $genotypesConditionalFile | tail -n +2) ) \
> phenotypesConditional$chrom.$region.txt 

phenoFile2=phenotypesConditional$chrom.$region.txt
head -n 3 $phenoFile2

# set outfile name
outFile=$basedir/data/GWAS/otherGWAS/$pheno/BOLTLMM.$vers/$pheno-BOLT-LMM-$vers-cond-chr${chrom}-reg$region-$start-$end.out

echo $outFile

# MEMORY NEEDED: 74,100,452,428 bytes = 5 C-NODES  <==== Sample QC SNPs
# MEMORY NEEDED: 96,739,226,940  = 7 C-NODES   <==== Oxford QC SNPs
# ALL CHROMOSOMES AT ONCE!!

# MUST MAKE SURE THAT THE .fam files for the autosomes and sex chroms are the same (they are as of May 16th)!!

# NOTE: This version uses 20 PCs as covariates - for LMM version as well!
# NOTE: This is the conditional version of $vers
# NOTE: consider using --covarUseMissingIndic option to include all samples, regardless of missing values in covariates. BOLT help:
#   --covarUseMissingIndic          include samples with missing covariates in 
#                                  analysis via missing indicator method 
#                                 (default: ignore such samples)

$bolt \
--bfile=$genoDataAuto \
--lmmInfOnly \
--modelSnps=$snpSet \
--exclude=snp.exclusion.$chrom.$region.txt \
--phenoFile=$phenoFile2 \
--phenoCol=$pheno \
--covarFile=$phenoFile2 \
--covarCol=Array \
--qCovarCol=Age.when.attended.assessment.centre \
--covarCol=Inferred.Gender \
--qCovarCol=PC{1:20} \
--qCovarCol=SNP{1:$nsnps} \
--numThreads=3 \
--statsFile=$outFile \
--covarMaxLevels=111 \
--geneticMapFile=/well/donnelly/ukbiobank_project_8874/clare/src/BOLT-LMM_v2.2/tables/genetic_map_hg19.txt.gz \
--verboseStats

rm $phenoFile2

# qsub -q short.qc -P donnelly.prjc runBOLT_versionNumber-LMM-conditional.sh

