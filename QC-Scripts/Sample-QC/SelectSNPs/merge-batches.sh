################################
# For UKBiobank sample QC pipeline step 1
################################
# Clare Bycroft, Jan 2016
#
# This script creates the base plink files required for Sample QC. Sex chromosomes and autosomes are created independently. There will be two sets: one with all batches combined, and another for each batch separately. All will contain the same set of QC'd SNPs (sex + autosome). This script merges all the batches together and also computes allele frequencies, which will probably be useful later.
################################

plink=/well/ukbiobank/qcoutput/Src/plink-v1.90b3k-x86_64/plink

############ %%%
inputdir=$basedir/data/ByBatch # this is the directory with plink files made by create-plink-subsets.sh
outputdir=$basedir/data/Combined
############



# SEX CHROMOSOMES

ls $inputdir/*sexchrom-sampleqc.bim | sed 's/\.bim//g' > plink-files-to-merge-sexchrom.txt
wc plink-files-to-merge-sexchrom.txt

$plink --merge-list plink-files-to-merge-sexchrom.txt --keep-allele-order --make-bed --out $outputdir/b1__b11-b001__b095-sexchrom-sampleqc

$plink --bfile $outputdir/b1__b11-b001__b095-sexchrom-sampleqc --keep-allele-order --freq --hardy --missing --out $outputdir/b1__b11-b001__b095-sexchrom-sampleqc



# AUTOSOMES

ls $inputdir/*autosome-sampleqc.bim | sed 's/\.bim//g' > plink-files-to-merge-autosome.txt
wc plink-files-to-merge-autosome.txt

$plink --merge-list plink-files-to-merge-autosome.txt --keep-allele-order --make-bed --out $outputdir/b1__b11-b001__b095-autosome-sampleqc

$plink --bfile $outputdir/b1__b11-b001__b095-autosome-sampleqc --keep-allele-order --freq --hardy --missing --out $outputdir/b1__b11-b001__b095-autosome-sampleqc



# PERFORM SOME CHECKS
echo "Were there any merge failures?" 
ls $outputdir | grep missnp

echo "N samples in *.fam files should be the same" 
wc $outputdir/b1__b11-b001__b095-autosome-sampleqc.fam
wc $outputdir/b1__b11-b001__b095-sexchrom-sampleqc.fam

diff $outputdir/b1__b11-b001__b095-autosome-sampleqc.fam $outputdir/b1__b11-b001__b095-sexchrom-sampleqc.fam 

echo "Do we still have all the SNPs?"
wc $outputdir/b1__b11-b001__b095-autosome-sampleqc.bim
wc $outputdir/b1__b11-b001__b095-sexchrom-sampleqc.bim
