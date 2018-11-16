#!/bin/bash
# The name of the job, can be anything, simply used when displaying the list of running jobs
#$ -N UKB-kinship-filtered
# Giving the name of the output log file
#$ -o Logs
#$ -j y
# One needs to tell the queue system to use the current directory as the working directory
# Or else the script may fail as it will execute in your top level home directory
#$ -cwd
#$ -V
#$ -P donnelly.prjc -q short.qc



plink=/well/ukbiobank/qcoutput/Src/plink-v1.90b3k-x86_64/plink
king14=/well/ukbiobank/qcoutput/Src/Linux-king-1.4/king
king19=/well/ukbiobank/qcoutput.V2_QCed.sample-QC/src/Linux-king-1.9/king

#NOTE: king19 is apparently faster, but will have to check against interim-release results before doing final run

## The input file is a list of the {nBatches choose 2} batches in the dataset:
## The first and the second columns specify the full path to the two batches,
## in plink binary format (all individuals, only SNPs for sample QC)
## The third column specifies the _directory_ to store the King results
basedir=/well/ukbiobank/qcoutput.V2_QCed.sample-QC

PairsOfBatches=b1__b11-b001__b095-pair_batches.txt
#PairsOfBatches=b001__b022-pair_batches.txt ## <----- for testing

lineno=$SGE_TASK_ID

Pair=`sed -n "$lineno{p;q}" $PairsOfBatches`
Batch1=$(echo $Pair | cut -f1 -d" ")
Batch2=$(echo $Pair | cut -f2 -d" ")
Merged=$(echo $Pair | cut -f3 -d" ")


# subset SNPs ---> filter for admixed samples 
############
threshold=0.003
############

SNPsToInclude=$basedir/QC-Scripts/PCA/pca-UKBio/b1__b11-b001__b095-autosome-sampleqc-fastpca-highquality-init-snpsToKeepPCs-$threshold.txt
# 93,512 SNPs

if [[ -f $Batch1.bed && -f $Batch1.bim && -f $Batch1.fam && -f $Batch2.bed && -f $Batch2.bim && -f $Batch2.fam ]] ; then
    
    if [[ ! -d $Merged ]] ; then
        mkdir $Merged
    fi
    
    ## Merge the two batches into a temporary plink binary dataset, i.e., generate $Merged.bed, $Merged.bim, $Merged.fam
    echo $plink --bfile $Batch1 --bmerge $Batch2.bed $Batch2.bim $Batch2.fam --allow-no-sex --make-bed --keep-allele-order --extract $SNPsToInclude --out $Merged
    time $plink --bfile $Batch1 --bmerge $Batch2.bed $Batch2.bim $Batch2.fam --allow-no-sex --make-bed --keep-allele-order --extract $SNPsToInclude --out $Merged
    
    ## Compute kinship coefficients of third degree or higher/closer
    ## $Merged is the name of the output directory (as King output several files)
    ## The relevant file would be $Merged/degree3.kin0
    #echo $king19 -b $Merged.bed --prefix $Merged/degree3-king19 --kinship --related --degree 3
    #time $king19 -b $Merged.bed --prefix $Merged/degree3-king19 --kinship --related --degree 3

    ## NOTE: only relevant when testing differences between the two verisons of KING
    #echo $king14 -b $Merged.bed --prefix $Merged/degree3-king14 --kinship --related --degree 3
    #time $king14 -b $Merged.bed --prefix $Merged/degree3-king14 --kinship --related --degree 3
    
    echo $king14 -b $Merged.bed --prefix $Merged/degree3-filtered --kinship --related --degree 3
    time $king14 -b $Merged.bed --prefix $Merged/degree3-filtered --kinship --related --degree 3
    
    ## Delete temporary plink files
    rm $Merged.*
    echo Done.
	
fi;
