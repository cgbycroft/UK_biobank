#!/bin/bash
#$ -cwd
#$ -N bolt_lmm-quant-Age-Ychrom_v4
#$ -o Logs
#$ -j y
#$ -pe shmem 3
#$ -V

# Downloaded from https://data.broadinstitute.org/alkesgroup/BOLT-LMM/downloads/ on May 12th 2016
bolt=/well/donnelly/ukbiobank_project_8874/clare/src/BOLT-LMM_v2.2/bolt 
basedir=/well/ukbiobank/qcoutput.V2_QCed.sample-QC

vers=v4         # as for v3 but sex chroms
pheno=X.mean.l2r   # actually X chromosome intensity! THIS IS FOR TESTING VARIANTS ON THE X CHROMOSOME; it assumes that the calls are still good....

### this sub-run is:
# BOLT version: LMM default - let BOLT decide on effect sizes.

mkdir $basedir/data/GWAS/YXintensity/BOLTLMM.$vers

genoDataAll=$basedir/data/Combined/byChrom/V2_QCed.export.UKBiLEVEAX_b1-b11.Batch_b001_b095-sexchroms # this has be subset using otherGWAS/subset-sex-chroms.sh

chr=$SGE_TASK_ID

thisChromGenoData=$genoDataAll-$chr

phenoFile=$basedir/data/GWAS/YXintensity/YchromPhenotypesForBOLT-v3.txt

# phenoFile and snpSet was created using subsetPhenotypes.R
snpSet=$basedir/data/PCA/b1__b11-b001__b095-autosome-sampleqc-fastpca-init-highquality.prune.in  # just use SNPs from original PCA :) NOT ACTUALLY USED IN LREG!!

outFile=$basedir/data/GWAS/YXintensity/BOLTLMM.$vers/Ychrom-BOLT-LMM-quant-Age-$vers-$chr.out

# MEMORY NEEDED: 74,100,452,428 bytes = 5 C-NODES  <==== Sample QC SNPs
# MEMORY NEEDED: 96,739,226,940  = 7 C-NODES   <==== Oxford QC SNPs

# NOTE: will probably need to do PCs 1-20 instead of 1-8 eventually! Need to add these to the phenotype file...

$bolt \
--bim=$thisChromGenoData-altcode.bim \
--bed=$thisChromGenoData.bed \
--fam=$thisChromGenoData.fam \
--phenoFile=$phenoFile \
--phenoCol=$pheno \
--covarFile=$phenoFile \
--qCovarCol=PC{1:20} \
--covarCol=Array \
--qCovarCol=Age.when.attended.assessment.centre \
--numThreads=3 \
--statsFile=$outFile \
--covarMaxLevels=111 \
--verboseStats


# change chromosome number in output back to sex chroms!

mv $outFile $outFile.orig
cat <(head -n 1 $outFile.orig) <( awk -v c=$chr '{OFS="\t"} {$2=c}1' $outFile.orig | tail -n +2) > $outFile


# qsub -t 23-26 -q short.qb -P donnelly.prjb runBOLT_v4-LMM.sh

