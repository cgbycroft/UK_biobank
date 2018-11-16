############ SAMPLE QC PIPELINE V1
# Here I am checking that the scripts for the pipeline for Release2 actually work when I run them
# I am using copies of the scripts in release1 qcoutput/QC-Scripts, where related to sample QC. Output files go to dummy output directories in the same folder.
# I will also note here which stages I have tested/completed etc. It will be a bit of a scratch-pad.
# This version will use plink files from:   /well/ukbiobank/expt/V2_QCed.SNP-QC/data/
# These files will have Oxford QC applied in this version, but is unlikely to be the final version. Data from previous versions are found in the .git repository under commit: "preliminary run v0"

#NOTE: version of R needs to be...


##############################
# 0. Initiation
##############################

basedir=/well/ukbiobank/qcoutput.V2_QCed.sample-QC
cd $basedir

# create list of batches for when using by-batch scripts
ls /well/ukbiobank/expt/V2_QCed.SNP-QC/data/ | grep -E '^Batch|^UK' > QC-Scripts/batchList.txt
batchList=$basedir/QC-Scripts/batchList.txt

# set list of helper scripts (e.g batch2sample.R) ---> use Colin's versions in final run!!
## Edit helper script "auxFunctions.R" to have correct base directories
helperscripts="/well/ukbiobank/expt/V2_QCed.SNP-QC/src/V2_QCed.snpqc-tests.R /well/ukbiobank/expt/V2_QCed.SNP-QC/src/V2_QCed.bin2clusterplots.R $basedir/QC-Scripts/R/scripts/readPSperformance.R $basedir/QC-Scripts/R/scripts/auxFunctions.R"

# just convenience for testing R scripts
#args = c("/well/ukbiobank/expt/V2_QCed.SNP-QC/src/V2_QCed.snpqc-tests.R","/well/ukbiobank/expt/V2_QCed.SNP-QC/src/V2_QCed.bin2clusterplots.R","/well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/R/scripts/readPSperformance.R","/well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/R/scripts/auxFunctions.R")

# create list of reference samples (takes a while so do this early)
Rscript $basedir/QC-Scripts/R/scripts/getReferences.R $helperscripts -out $basedir/QC-Scripts/referenceList.txt

# create list of samples (bestarray) and their respective batch.
Rscript $basedir/QC-Scripts/R/scripts/list-sample-batches.R $helperscripts -out $basedir/QC-Scripts/sampleBatchList.txt

# get number of batches (for array jobs)
nbatches=`cut -f2 -d' ' $basedir/QC-Scripts/sampleBatchList.txt | sort -u | wc -l | cut -f1 -d' '`



##############################
# 1. Create plink files for sample QC; compute basic stats
##############################

cd QC-Scripts/Sample-QC/SelectSNPs

# select snps on autosomes and sex chroms.
Rscript selectSNPS.R $helperscripts > Logs/selectSNPS.log &

# run some plots on snp selection results (%%% not complete, but not required for pipeline)
Rscript plot-selectSNPs.R $helperscripts > Logs/plot-selectSNPS.log & 

# apply subsets by batch (also creates plink files with reference samples)
qsub -t 1-$nbatches ./create-plink-subsets.sh

# merge batches (after above has run)
source ./merge-batches.sh > Logs/merge-batches.log &

# Compute basic statistics
source ./compute-basic-stats.sh > Logs/compute-basic-stats.log &

# APPLY OXFORD QC PER BATCH (not for Sample QC, but for GWAS tests ===> but ultimately I used the data that Colin created.)

# Create combined file with SNP exclusions by batch 
qsub -t 1-$nbatches ./create-plink-subsets-Oxford-QC.sh

# merge batches (after above has run) and split by chromosome (for GWAS)
source ./merge-batches-oxfordqc.sh > Logs/merge-batches-oxfordqc.log &



##############################
# 2. Relatedness
##############################

cd $basedir/QC-Scripts/Sample-QC/Relatedness

# Create the file with the list of batch-pairs.  %%% Will need to update this file in the final run
source ./list-unique-pairs.sh

## Run King in pairs. First round, just use all the SNPs.
# make sure reference samples are removed from plink base files (otherwise a waste of computing power)
npairs=`wc -l b1__b11-b001__b095-pair_batches.txt | cut -f1 -d' '`
qsub -t 1-$npairs submit-pt1-kinship-in-pairs.sh

# combine king files across batches.
source ./submit-pt2-combine-kin-files.sh

# plot initial round results
Rscript plot-king.R $helperscripts -in $basedir/data/Relatedness/b1__b11-b001__b095-pair_batches-init.kin0 -ibs0 0.002 -out plots > Logs/plot-king-init.log 

# prune the het-miss outliers
Rscript postking/filter-king.R $helperscripts -in $basedir/data/Relatedness/b1__b11-b001__b095-pair_batches-init.kin0 -ibs0 0.002 -outliers $basedir/QC-Scripts/Sample-QC/HetMissing/b1__b11-b001__b095-autosome-sampleqc-hetcorrected-6PCs-imiss-all-hetmiss-outliers.txt > Logs/filter-king-init.log &

# find maximal set of unrelated samples; (NOTE: set seed, as there is a random element to this algorithm)
Rscript postking/get-unrelated.R $helperscripts -in $basedir/data/Relatedness/b1__b11-b001__b095-pair_batches-init.kin0 -ibs0 0.002 -outdir postking > Logs/get-unrelated-init.log 

# how many unrelated samples?
grep Total Logs/get-unrelated-init.log


# ====> NOW RUN PCA


cd $basedir/QC-Scripts/Sample-QC/Relatedness

# KING filtered round
qsub -t 1-$npairs submit-pt1-kinship-in-pairs-filtered.sh

# combine files
source ./submit-pt2-combine-kin-files-filtered.sh > Logs/submit-pt2-combine-kin-files-filtered.log

# plot results
Rscript plot-king.R $helperscripts -in $basedir/data/Relatedness/b1__b11-b001__b095-pair_batches.filtered.kin0 -ibs0 0.0012 -out plots > Logs/plot-king-filtered.log


# ====> NOW RUN HET/MISSING 


Rscript postking/filter-king.R $helperscripts -in $basedir/data/Relatedness/b1__b11-b001__b095-pair_batches.filtered.kin0 -ibs0 0.0012 -outliers $basedir/QC-Scripts/Sample-QC/HetMissing/b1__b11-b001__b095-autosome-sampleqc-hetcorrected-6PCs-imiss-all-hetmiss-outliers.txt > Logs/filter-king.log &

# plot the fully pruned set (this is pruned to exclude >10 3rd-degree relatives too!)
Rscript plot-king.R $helperscripts -in $basedir/data/Relatedness/b1__b11-b001__b095-pair_batches.filtered-samples-pruned.kin0 -ibs0 0.0012 -out plots > Logs/plot-king-filtered-pruned.log

# plot the middle pruned set (this is pruned to exclude het-miss outliers + anyone with >200 relatives
Rscript plot-king.R $helperscripts -in $basedir/data/Relatedness/b1__b11-b001__b095-pair_batches.filtered-samples-pruned-200.kin0 -ibs0 0.0012 -out plots > Logs/plot-king-filtered-pruned-200.log

# plot the least pruned set -- only excluding het-miss outliers
Rscript plot-king.R $helperscripts -in $basedir/data/Relatedness/b1__b11-b001__b095-pair_batches.filtered-hetmiss-out.kin0 -ibs0 0.0012 -out plots > Logs/plot-king-filtered-hetmiss-out.log

# find families (make sure to first remove *families.RData if the above has been re-done!)
Rscript postking/find-families.R $helperscripts -in $basedir/data/Relatedness/b1__b11-b001__b095-pair_batches.filtered-samples-pruned.kin0 -ibs0 0.0012 > Logs/find-families-filtered-pruned.log

# find families fully pruned set excluding duplicates
Rscript postking/find-families.R $helperscripts -in $basedir/data/Relatedness/b1__b11-b001__b095-pair_batches.filtered-samples-pruned.kin0 -ibs0 0.0012 -exclude-dupes > Logs/find-families-filtered-pruned-nodupes.log

# find families middle pruned set (keep duplicates) <====== THIS WAS THE SET WE USED TO FIND DUPLICATES
Rscript postking/find-families.R $helperscripts -in $basedir/data/Relatedness/b1__b11-b001__b095-pair_batches.filtered-samples-pruned-200.kin0 -ibs0 0.0012 > Logs/find-families-filtered-pruned-200.log

# find families middle pruned set (exclude duplicates) <====== THIS WAS THE SET WE RELEASED
Rscript postking/find-families.R $helperscripts -in $basedir/data/Relatedness/b1__b11-b001__b095-pair_batches.filtered-samples-pruned-200.kin0 -exclude-dupes -ibs0 0.0012 > Logs/find-families-filtered-pruned-200-nodupes.log

# find families in actually released kinship table (excludes duplicates and some other ad-hocs)
Rscript postking/find-families.R $helperscripts -in $basedir/data/ForRelease/b1__b11-b001__b095-sampleTable_v4_Kinship.txt -ibs0 0.0012 > Logs/find-families-sampleTable_v4_Kinship.log

# find families in actually released kinship table (excludes duplicates and some other ad-hocs) + exclude kinship < 0.05
Rscript postking/find-families.R $helperscripts -in $basedir/data/ForRelease/b1__b11-b001__b095-sampleTable_v3_Kinship.txt -ibs0 0.0012 -3degree 0.05 > Logs/find-families-sampleTable_v3_Kinship-3deg-0.05.log

# plot kinship in actually released kinship table (excludes duplicates and some other ad-hocs)
Rscript plot-king.R $helperscripts -in $basedir/data/ForRelease/b1__b11-b001__b095-sampleTable_v4_Kinship.txt -ibs0 0.0012 -out plots > Logs/plot-king-sampleTable_v4_Kinship.log &


# Find trio lists ====> For Jonathan.
Rscript postking/find-trios.R $helperscripts -in $basedir/data/Relatedness/b1__b11-b001__b095-pair_batches.filtered-samples-pruned.kin0 -ibs0 0.0012 -exclude-dupes > Logs/find-trios-filtered-pruned-nodupes.log

Rscript postking/find-trios.R $helperscripts -in $basedir/data/ForRelease/b1__b11-b001__b095-sampleTable_v4_Kinship.txt -ibs0 0.0012 > Logs/find-trios-sampleTable_v4_Kinship.log


# Find maximal set of unrelated samples (look at plots to set ibs0 threshold) (NOTE: set seed???, as there is a random element to this algorithm)
Rscript postking/get-unrelated.R $helperscripts -in $basedir/data/Relatedness/b1__b11-b001__b095-pair_batches.filtered.kin0 -ibs0 0.0012 \
-exclude $basedir/QC-Scripts/Sample-QC/HetMissing/b1__b11-b001__b095-autosome-sampleqc-hetcorrected-6PCs-imiss-all-hetmiss-outliers.txt \
-exclude $basedir/QC-Scripts/Sample-QC/Relatedness/b1__b11-b001__b095-pair_batches.filtered-excess-relatives.txt \
-outdir postking > Logs/get-unrelated-filtered.log


# =====> NOW RUN PCA ROUND 2, using the list of unrelated and high quality samples made above.

####################
# 2.b COMPARISON WITH PLINK IBD (not fundamental to pipeline)
####################

# ====> After selecting WhiteBritish set, and pruning kinship table, run plink IBD to validate pruned kinship estimates.
Rscript select-IBD-samples.R $helperscripts -wb $basedir/QC-Scripts/WhiteBritish/b1__b11-b001__b095-pca-UKbio-round2-White_British.txt -rel $basedir/data/Relatedness/b1__b11-b001__b095-pair_batches.filtered-samples-pruned-200.kin0 > Logs/select-IBD-samples.log

qsub submit-plink-ibd.sh   # parallel option (works - takes about 2hrs)

# create combined file with unique pairs, such that it looks like the KING output.
Rscript process-IBD-Kinship.R $helperscripts -ibd $basedir/data/Relatedness/b1__b11-b001__b095-autosome-sampleqc-filtered.0.003-White_British.genome -hm $basedir/QC-Scripts/Sample-QC/HetMissing/b1__b11-b001__b095-autosome-sampleqc-hetcorrected-6PCs-imiss-all-hetmiss-outliers.txt > Logs/process-IBD-Kinship.log

# filter plink ibs like the other king results
Rscript postking/filter-king.R $helperscripts -in $basedir/data/Relatedness/b1__b11-b001__b095-autosome-sampleqc-filtered.0.003-White_British.genome.combined -ibs0 0.0012 -outliers $basedir/QC-Scripts/Sample-QC/HetMissing/b1__b11-b001__b095-autosome-sampleqc-hetcorrected-6PCs-imiss-all-hetmiss-outliers.txt > Logs/filter-plink-ibd.log

# plink results (no pruning)
Rscript plot-king.R $helperscripts -in $basedir/data/Relatedness/b1__b11-b001__b095-autosome-sampleqc-filtered.0.003-White_British.genome.combined -ibs0 0.0012 -out plots > Logs/plot-plink-ibd.log

# plink results (hetmiss exclusions only)
Rscript plot-king.R $helperscripts -in $basedir/data/Relatedness/b1__b11-b001__b095-autosome-sampleqc-filtered.0.003-White_British.genome.combined-hetmiss-out.kin0 -ibs0 0.0012 -out plots > Logs/plot-plink-ibd-hetmiss-out.log

# plink results (white british only + removing > 200 relatives)
Rscript plot-king.R $helperscripts -in $basedir/data/Relatedness/b1__b11-b001__b095-autosome-sampleqc-filtered.0.003-White_British.genome.combined-samples-pruned-200.kin0 -ibs0 0.0012 -out plots > Logs/plot-plink-pruned-ibd-200.log

# King results (white british only)
Rscript plot-king.R $helperscripts -in $basedir/data/Relatedness/b1__b11-b001__b095-pair_batches.filtered-samples-pruned-200-White_British.kin0 -ibs0 0.0012 -out plots > Logs/plot-king-filtered-pruned-200-White_British.log


# Compare King and Plink (for Paper only)
Rscript compare-king-plink.R $helperscripts -ibd $basedir/data/Relatedness/b1__b11-b001__b095-autosome-sampleqc-filtered.0.003-White_British.genome.combined-samples-pruned-200.kin0 -rel /well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/Relatedness/b1__b11-b001__b095-pair_batches.filtered-samples-pruned-200-White_British.kin0 -ibs0 0.0012 > Logs/compare-king-plink1.log


####################
# 2.c CONCORDANCE BETWEEN BLIND-SPIKED DUPLICATES
####################

cd $basedir/QC-Scripts/Sample-QC/Relatedness

nDupes=`tail -n +2 b1__b11-b001__b095-pair_batches.filtered-samples-pruned-200-duplicates-twins.txt | wc -l | cut -f1 -d' '`

# Using SNPs only used in sample QC
qsub -t 1-$nDupes -P donnelly.prja -q short.qa ./postking/duplicate-concordance.sh

# Using all SNPs left after applying our QC (just run all duplicates together on the server. It is faster.)
source ./postking/duplicate-concordance-release.sh

#rm keep.*.txt


# run based on pruned kinship table
Rscript ./postking/plot-duplicate-concordance.R $helperscripts -in b1__b11-b001__b095-pair_batches.filtered-samples-pruned-200 > Logs/plot-duplicate-concordance.filtered-samples-pruned-200.log &

# run based on unpruned kinship table
Rscript ./postking/plot-duplicate-concordance.R $helperscripts -in b1__b11-b001__b095-pair_batches.filtered > Logs/plot-duplicate-concordance.filtered.log &

# Compute the differences and save lists of discordant snps (added Aug 2018)
K=`wc -l b1__b11-b001__b095-pair_batches.filtered-duplicates-twins.columnsInGenotypeFile.txt | cut -f1 -d' '`
qsub -t 1-$K -P donnelly.prjc -q short.qc -N dup-concordance-bySNP ./postking/submit-duplicate-bySNPdifferences.sh



##############################
# 3. PCA
##############################

cd $basedir/QC-Scripts/PCA/pca-UKBio

# create list of samples to exclude, using initial round of KING
Rscript pca-sample-filters.R $helperscripts -k $basedir/data/Relatedness/b1__b11-b001__b095-pair_batches-init.kin0 -ibs0 0.002 -out pca-UKBio-sample-exclusions-init.txt  > Logs/pca-sample-filters-init.log

# subset the data by samples and SNPs; prune for LD
source ./submit-pca-subset-init.sh > Logs/submit-pca-subset-init.log

# submit job for computing PCA
qsub -q long.qc -P donnelly.prjc submit-fastpca-init.sh

# run projections 
qsub -t 1-$nbatches submit-pca-UKBio-init.sh 

# plot PCs and create RData file
Rscript plot-pca-UKBio.R $helperscripts -in $basedir/data/PCA/b1__b11-b001__b095-autosome-sampleqc-fastpca-highquality-init.snpload.map -out b1__b11-b001__b095-pca-UKbio-init

# filter SNPs for KING
Rscript filter-snps-for-king.R $helperscripts -in b1__b11-b001__b095-autosome-sampleqc-fastpca-highquality-init -threshold 0.003 >  Logs/filter-snps-for-king-init.log

# ====> Now run KING for the second round + HET/MISSING

# PCA round 2
# create list of samples to exclude, using second round of KING and het-missing list
Rscript pca-sample-filters.R $helperscripts -k $basedir/data/Relatedness/b1__b11-b001__b095-pair_batches.filtered.kin0 -ibs0 0.0012 \
-exclude $basedir/QC-Scripts/Sample-QC/HetMissing/b1__b11-b001__b095-autosome-sampleqc-hetcorrected-6PCs-imiss-all-hetmiss-outliers.txt \
-exclude $basedir/QC-Scripts/Sample-QC/Relatedness/b1__b11-b001__b095-pair_batches.filtered-excess-relatives.txt \
-out pca-UKBio-sample-exclusions-round2.txt  > Logs/pca-sample-filters-round2.log

# subset the data by samples and SNPs; prune for LD
source ./submit-pca-subset-round2.sh > Logs/submit-pca-subset-round2.log

# submit job for computing PCA
#qsub -q short.qc -P donnelly.prjc submit-fastpca-round2.sh
qsub -q coolibah.q -P donnelly.prjc submit-fastpca-round2-coolibah.sh

# run projections 
qsub -t 1-$nbatches submit-pca-UKBio-round2.sh

# plot the results (note: this must be run immediately after projections are created, as projection file names are not unique)
Rscript plot-pca-UKBio.R $helperscripts -in $basedir/data/PCA/b1__b11-b001__b095-autosome-sampleqc-fastpca-highquality-round2-coolibah.snpload.map -out b1__b11-b001__b095-pca-UKbio-round2

# ====> Now select White-British subset(s)

# PCA on white-british subset <======= THIS IS WHAT I USED IN THE GWASs
Rscript pca-sample-filters.R $helperscripts -k $basedir/data/Relatedness/b1__b11-b001__b095-pair_batches.filtered.kin0 -ibs0 0.0012 \
-exclude $basedir/QC-Scripts/Sample-QC/HetMissing/b1__b11-b001__b095-autosome-sampleqc-hetcorrected-6PCs-imiss-all-hetmiss-outliers.txt \
-exclude $basedir/QC-Scripts/Sample-QC/Relatedness/b1__b11-b001__b095-pair_batches.filtered-excess-relatives.txt \
-exclude $basedir/QC-Scripts/WhiteBritish/b1__b11-b001__b095-pca-UKbio-round2-White_British-exclusions.txt \
-out pca-UKBio-sample-exclusions-White_British.txt > Logs/pca-sample-filters-White_British.log

source ./submit-pca-subset-White_British.sh > Logs/submit-pca-subset-White_British.log 
qsub -q coolibah.q -P donnelly.prjc submit-fastpca-White_British.sh

# project the rest of the samples (although this doesn't really make sense as we only used WB inds in the PCA)
qsub -t 1-$nbatches submit-pca-UKBio-White_British.sh
Rscript plot-pca-UKBio.R $helperscripts -in $basedir/data/PCA/b1__b11-b001__b095-autosome-sampleqc-fastpca-highquality-White_British.snpload.map -out b1__b11-b001__b095-pca-UKbio-White_British

# just plot the subset of WB samples
Rscript plot-pca-UKBio.R $helperscripts -in $basedir/data/PCA/b1__b11-b001__b095-autosome-sampleqc-fastpca-highquality-White_British.snpload.map -out b1__b11-b001__b095-pca-UKbio-White_British_subset


############
# White British PCA using more rare variants <== for thesis
source ./submit-pca-subset-White_British_rare.sh > Logs/submit-pca-subset-White_British_rare.log &
qsub -q coolibah.q -P donnelly.prjc submit-fastpca-White_British_rare.sh

# project the rest of the samples (although this doesn't really make sense as we only used WB inds in the PCA)
qsub -t 1-$nbatches submit-pca-UKBio-White_British_rare.sh
# TO HERE- 24/8/2017

Rscript plot-pca-UKBio.R $helperscripts -in $basedir/data/PCA/b1__b11-b001__b095-autosome-sampleqc-fastpca-highquality-White_British_rare.snpload.map -out b1__b11-b001__b095-pca-UKbio-White_British_rare


# just plot the subset of WB samples (and save the PCs as an R file!!)
Rscript plot-pca-UKBio.R $helperscripts -in $basedir/data/PCA/b1__b11-b001__b095-autosome-sampleqc-fastpca-highquality-White_British_rare.snpload.map -out b1__b11-b001__b095-pca-UKbio-White_British_subset_rare
# More plotting (e.g on a map is done in plot-White_British-PCA_v2.R)





##############################
# 4. Het/missing
##############################
# NOTE: missing rates and heterozygosity are computed first on all snps and samples when creating plink files, as these feed into the PCs. Then after computing PCs one can adjust the heterozygosity for ancestry - this is what is done here. The script compute-het-missing.sh only computes it over the set of SNPs used in the final Kinship calculation. Otherwise we use the output from compute-basic-stats.sh

cd $basedir/QC-Scripts/Sample-QC/HetMissing

# ~100 mins each
# plink het and missing computation over combined plink base files. %%% check input files are correct.
source ./compute-het-missing.sh > Logs/compute-het-missing.log &

# LROH
source ./compute-LROH.sh > Logs/compute-LROH.sh & # <====== IT RUNS

# plotting ( only use this if you don't already have the PCs to do the corrections - i.e if you haven't run plot-pca-UKBio.R )
Rscript plot-het-missing-uncorrected.R $helperscripts -in $basedir/data/Combined/b1__b11-b001__b095-autosome-sampleqc

nPCsadjust=6PCs

# Adjust with PCs
Rscript compute-het-corrected.R $helperscripts -in $basedir/data/Combined/b1__b11-b001__b095-autosome-sampleqc -pcs $basedir/data/PCA/b1__b11-b001__b095-pca-UKbio-init.RData -npcs $nPCsadjust

# plot adjusted and unadjusted
Rscript plot-het-missing-corrected.R $helperscripts -in $basedir/data/Combined/b1__b11-b001__b095-autosome-sampleqc -npcs $nPCsadjust


# run british aberrant test (uses markdown language)
#Rscript -e "require(knitr); knit2html('hetmiss-British-analysis.Rmd')" $basedir/data/Combined/b1__b11-b001__b095-autosome-sampleqc-hetcorrected-$nPCsadjust-imiss

# run british + Irish + Indian + Other white - aberrant (uses markdown language). Make sure to change EVAL if you want to re-run the aberrant analysis.
Rscript -e "require(knitr); knit2html('hetmiss-British-Irish-White-Indian-analysis.Rmd')" $basedir/data/Combined/b1__b11-b001__b095-autosome-sampleqc-hetcorrected-$nPCsadjust-imiss

# add missingness criteria for non BIWI samples and compile a complete list of outliers
Rscript compile-het-missing-exclusion-list.R $helperscripts -in $basedir/data/Combined/b1__b11-b001__b095-autosome-sampleqc-hetcorrected-$nPCsadjust-imiss > Logs/compile-het-missing-exclusion-list.log

# =====> NOW RUN KINSHIP FILTERING, AND PCA ROUND 2

##############################
# 5. Select ‘Caucasian’ subset
##############################
# NOTE: This requires PCs to have been computed for all samples

cd $basedir/QC-Scripts/WhiteBritish
# script to create BatchesFile with list of file names for the projected pcs in each batch (or just incorporate it into the subset-White-British.R code??

# run aberrant (make sure the right R version is used)
Rscript subset-White-British.R $helperscripts -in b1__b11-b001__b095-pca-UKbio-round2 > Logs/subset-White-British.log & # <====== 

# =====> Then run PCA on just these samples. Choose a set of SNPs again (so to capture variation that is rare in all samples but not in this set - probably mostly the same SNPs, but might pay to check). Also remove related samples within this subset.


##############################
# 6. Sex mismatches/aneuploidies
##############################

### 6a. Extract and summarise CNV data for the sex chromosomes
cd $basedir/QC-Scripts/Sample-QC/ExtractCNV

# get bytes lengths for each batch file ===> didn't actually use this!!
#qsub -t 1-$nbatches -P donnelly.prja -q short.qa submit-file-indexes.sh

# get order of SNPs in cnv file (takes about 5mins each)
awk '{print $1,$2,$3}' <(zcat $basedir/../source/V2_All/Batch_b001/AxiomGT1.cnv.log2ratio.final.txt.gz) > Batch_b001-CNV-snp-order.txt
awk '{print $1,$2,$3}' <(zcat $basedir/../source/V2_All/UKBiLEVEAX_b1/AxiomGT1.cnv.log2ratio.final.txt.gz) > UKBiLEVEAX_b1-CNV-snp-order.txt

# submit scripts for by-batch extraction of sex chromosome intensities
qsub -t 12-$nbatches -P donnelly.prjc -q short.qc submit-extract-sex-CNV.sh

# convert text files to hdf5 (takes a while...)
Rscript convert-to-hdf5.R $helperscripts

# compute mean & sd values for baf and l2ratio across each sex chromosome
Rscript compute-cnv-summaries.R $helperscripts -out b1__b11-b001__b095-sexchrom-sampleqc-cnv-summaries

### 6b. look at various properties of sex chromosomes to decide in list to exclude for phasing
cd $basedir/QC-Scripts/Sample-QC/sexChroms

# look at properties of snps across sex chromosomes
Rscript plot-snpqc-sexchroms.R $helperscripts

# Adjust X-heterozygosity with PCs
nPCsadjust=6PCs
Rscript ../HetMissing/compute-het-corrected.R $helperscripts -in $basedir/data/Combined/b1__b11-b001__b095-sexchrom-sampleqc-23 -pcs $basedir/data/PCA/b1__b11-b001__b095-pca-UKbio-init.RData -npcs $nPCsadjust

# plot adjusted and unadjusted
Rscript plot-het-missing-corrected.R $helperscripts -in $basedir/data/Combined/b1__b11-b001__b095-autosome-sampleqc -npcs $nPCsadjust

# find list of samples with odd sex information (and plot various statistics)
Rscript gender-checks.R $helperscripts -out b1__b11-b001__b095-sexchrom-sampleqc-sexCheck




##############################
# 7. Create lists of exclusion with flags for researchers + phasing
##############################

cd $basedir/QC-Scripts/Sample-QC/Flags

# NOTE: duplicate list is from Colin, after considering Sam Murphy's list of twins checked by UKBiobank.
Rscript create-exclusion-flags.R $helperscripts \
-hm $basedir/QC-Scripts/Sample-QC/HetMissing/b1__b11-b001__b095-autosome-sampleqc-hetcorrected-6PCs-imiss-all-hetmiss-outliers.txt \
-rel1 $basedir/data/Relatedness/b1__b11-b001__b095-pair_batches.filtered-samples-pruned-200-sample-gt10_3rdDeg_relatives.txt \
-rel2 $basedir/data/Relatedness/b1__b11-b001__b095-pair_batches.filtered-samples-pruned-200-sample-effective_exclusions.txt \
-sex $basedir/QC-Scripts/Sample-QC/sexChroms/b1__b11-b001__b095-sexchrom-sampleqc-sexCheck-phaseExclusionsCriteria1.txt \
-dup /well/ukbiobank/expt/V2_QCed.identical_samples/data/V2_QCed.duplicates_exclude.txt \
> Logs/create-exclusion-flags.log



##############################
# 8. Create sample table for researchers (takes in flags file)
##############################

mkdir $basedir/data/ForRelease

Rscript create-sample-table.R $helperscripts \
-qcflags $basedir/QC-Scripts/Sample-QC/Flags/b1__b11-b001__b095-sampleQC-flags.txt \
-pcs /well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/PCA/b1__b11-b001__b095-pca-UKbio-round2 \
-wb $basedir/QC-Scripts/WhiteBritish/b1__b11-b001__b095-pca-UKbio-round2-White_British.txt \
-fam /well/ukbiobank/expt/V2_QCed.export/data/calls_export/v3/V2_QCed.export.UKBiLEVEAX_b1-b11.Batch_b001_b095.chr22.fam \
-het $basedir/data/Combined/b1__b11-b001__b095-autosome-sampleqc-hetcorrected-imiss.RData \
-rel $basedir/data/Relatedness/b1__b11-b001__b095-pair_batches.filtered-samples-pruned-200.kin0 \
-out b1__b11-b001__b095-sampleTable_v4 \
> Logs/create-sample-table_v4.log &

# v3: all released columns (Dec 2016) ==> sent to Alan
# v4: as v3 but excluding 33 withdrawn samples (Feb 2017). Uses v3 of genotype calls.



##############################
# FROM HERE ON IS ANALYSIS THAT IS NOT GENERATING OUTPUT FOR RESEARCHERS
##############################


##############################
# 9. Compare genotype calls to interim release
##############################

cd $basedir/QC-Scripts/Sample-QC/InterimComparison

source ./compute-frequencies.sh # compute frequencies based on final data.

qsub ./process-plink-files.sh
qsub ./count-differences2.sh
    
Rscript plot-comparison.R # run in live session


##############################
# 10. Compare genotype calls to ExAC where they overlap
##############################

cd $basedir/QC-Scripts/ExacComparison

source ./extract-AF.sh > Logs/extract-AF.log

Rscript process-Exac-UKBB.R $helperscripts -exac $basedir/data/Exac/ExAC.r1.sites.vep.UKBB-SNPs.aCounts_NFE.txt -samples b1__b11-b001__b095-sampleTable_v4




##############################
# 11. Get numbers for tables and/or figures for paper
##############################

cd $basedir/forPaper

# Just run in a live session
Rscript $basedir/QC-Scripts/forPaper/numbersForPaper.R -in b1__b11-b001__b095-sampleTable_v4



##############################
# NA. cluster plots
##############################



##############################
# NA. Git repository (run approx. every week)
##############################
cd /well/ukbiobank/qcoutput.V2_QCed.sample-QC

git add ./*.sh  # <--- THESE COMMANDS DON"T ACTUALLY FIND EVERYTHING!!!!
git add ./*.R
git add ./*.txt
#git add ./\*.png
#git add data/Relatedness/b1__b11-b001__b095-pair_batches.filtered.kin0
#git add data/Relatedness/b1__b11-b001__b095-pair_batches-init.kin0
#git add data/PCA/b1__b11-b001__b095-pca-UKbio-init.RData

find . -name '*.R' | xargs git add  # <---- INSTEAD USE THESE AND CHECK checkStatus.tmp
find . -name '*.sh' | xargs git add
#find . -name '*.txt' | xargs git add
find . -name '*.py' | xargs git add

git status > checkStatus.tmp

#git status

#git commit -m "KING and PCA initial runs. Create QCd set for GWAS. Version 1 of GWAS tests"
#git commit -m "Updated KING plotting functions; added new het-missing scripts and run them. Some updates too GWAS tests."
#git commit -m "Mainly scripts to check/plot output from GWAS testing. v3 of GWAS."
#git commit -m "Mainly scripts for sex chromosome checks and v5 of GWAS."
#git commit -m "Mainly scripts for sex chromosome HMM and v6 of GWAS. Also GWASpipelineTemplate created."
#git commit -m "Updates for sex chromosome HMM."
#git commit -m "new Y chroms GWAS + filtering for X in phasing"
#git commit -m "Re-ran kinship filtering (same output as before + extra files with other lists). Updated create-sample-table.R file and re-ran these to create release data."
#git commit -m "Re-ran kinship file creation for release. Added kinship field descriptions file. Created scripts for plink ibd calculation and plotted plink ibd white british only compared to King white british only. Created new sample table (v3) with extra kinship flag."
#git commit -m "Wrote and ran find-families.R and find-trios.R on output kinship table. Wrote and ran script for duplicates concordance. Created new sample table (v4) excluding 33 withdrawl samples."
#git commit -m "Updates of most scripts. New folder for plots and numbers for paper. "
#git commit -m "Data and script freeze. Includes some text files."
#git commit -m "Updated scripts for figures and numbers for paper; exac analysis and interim comparison; new version of cluster plotting scripts to read in data for final release."
#git commit -m "X chromosome GWAS."
#git commit -m "Extra checks and plots for X chromosome. Sex changes from interim check for Alan."
#git commit -m "Changes to plotting for paper. Plus the credible set analysis for paper."
git commit -m "Adding readme for repository online."


