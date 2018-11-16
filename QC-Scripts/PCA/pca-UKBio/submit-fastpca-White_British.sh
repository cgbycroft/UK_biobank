#!/bin/bash
# The name of the job, can be anything, simply used when displaying the list of running jobs
#$ -N UKBB-pca-UKBio-fastpca-White_British
# Giving the name of the output log file
#$ -o Logs
#$ -j y
# One needs to tell the queue system to use the current directory as the working directory
# Or else the script may fail as it will execute in your top level home directory
#$ -cwd
#$ -V
## Ask for 8 cores = 128GB
#$ -pe shmem 8


## smartpca requires libgsl.so.0
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/gpfs0/apps/well/gsl/1.16-gcc4.7.2/lib


smartpca=/well/ukbiobank/qcoutput/Src/EIG6.0.1/bin/smartpca
snploads=/well/ukbiobank/qcoutput/Src/coolibah/snploads/snploads
nPCs=40

basedir=/well/ukbiobank/qcoutput.V2_QCed.sample-QC

#######################
## I have picked one batch in order to test the script
## In general, use plink to merge many (all?) batches,
## at the autosomal SNPs only.
##
## I will also assume that SNP+sample QC has already been applied,
## that SNPs have been pruned for high LD and the specific regions
## of very high LD have been removed, as described in 
## $basedir/QC-Scripts/pca-UKBio/flashpca-example/submit-flashpca.sh
Input=$basedir/data/PCA/b1__b11-b001__b095-autosome-sampleqc-fastpca-White_British-highquality-pruned
Output=$basedir/data/PCA/b1__b11-b001__b095-autosome-sampleqc-fastpca-highquality-White_British
#######################


cp $Input.bim $Input.pedsnp  # <========== add this back in!

cat -n $Input.fam | awk '{print $1,$1,$4,$5,$6,$7}' > $Input.pedind # <========== add this back in!

## Actually smartpca takes just one argument and it is the name
## of a parameter file, which I have called "par.smartpca.fast"

printf "genotypename: $Input.bed
snpname: $Input.pedsnp
indivname: $Input.pedind
evecoutname: $Output.evecs
fastmode: YES
numoutevec: $nPCs
fastdim: 50
fastiter: $nPCs
" > par.smartpca.White_British.fast
# other possible parameters include 
#fastdim:   (default 2 * numeigs) number of dimensions in fast pca approximation
#fastiter:  (default  numeigs)    number of iterations in fast pca approximation
#numthreads:   10 # although apparently this doesn't apply for fastmode=YES...

time $smartpca -p par.smartpca.White_British.fast



##############################################                                                                                                                                
## FastPCA is more efficient than FlashPCA but does not compute SNP loadings                                                                                                  
## Here is how to do it with a small program that uses Eigen                                                                                                                  
## Remove the first line, the first column and the last column,                                                                                                               
## to leave just the principal components as a large matrix of numbers only                                                                                                
awk '{$(NF--)=""; $1=""; print}' $Output.evecs | sed '1d' > $Output.evecs.tempi


## It is important that this plink dataset is the same dataset that FastPCA used                                                                                              
## to compute the principal components in the first place                                                                                                                     
time $snploads --bfile $Input --projpath $Output.evecs.tempi --snploads $Output.snpload


#rm $Output.evecs.tempi

# This next bit is the same as in the flashpca script
## Create shellfish-style *.snpload.map file for the precomputed loadings
## It is crucial to use the keep-allele-order option so that MAF is the frequency of Allele 1
## If we want to project new data, we have to standardize the genetic data in exactly the
## same as it was standardized for the PCA computation
## You can check whether new data is projected correctly by projecting the samples used in 
## the principal component analysis -- are the projections equal to the principal components?

$plink --bfile $Input --freq --keep-allele-order --out $Output
awk '{print $5}' $Output.frq | awk 'NR>1' > $Output.mean
paste $Input.bim $Output.mean $Output.snpload > $Output.snpload.map


## run this on short qc. qsub -q coolibah.q -P donnelly.prjc submit-fastpca-round2.sh
