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

ls $inputdir/*sexchrom-oxfordqc.bim | sed 's/\.bim//g' > plink-files-to-merge-oxfordqc-sexchrom.txt
wc plink-files-to-merge-oxfordqc-sexchrom.txt

$plink --merge-list plink-files-to-merge-oxfordqc-sexchrom.txt --keep-allele-order --make-bed --out $outputdir/b1__b11-b001__b095-sexchrom-oxfordqc

$plink --bfile $outputdir/b1__b11-b001__b095-sexchrom-oxfordqc --keep-allele-order --freq --hardy --missing --out $outputdir/b1__b11-b001__b095-sexchrom-oxfordqc



# AUTOSOMES

ls $inputdir/*autosome-oxfordqc.bim | sed 's/\.bim//g' > plink-files-to-merge-oxfordqc-autosome.txt
wc plink-files-to-merge-oxfordqc-autosome.txt

$plink --merge-list plink-files-to-merge-oxfordqc-autosome.txt --keep-allele-order --make-bed --out $outputdir/b1__b11-b001__b095-autosome-oxfordqc

$plink --bfile $outputdir/b1__b11-b001__b095-autosome-oxfordqc --keep-allele-order --freq --hardy --missing --out $outputdir/b1__b11-b001__b095-autosome-oxfordqc



# PERFORM SOME CHECKS
echo "Were there any merge failures?" 
ls $outputdir | grep missnp

echo "N samples in *.fam files should be the same" 
wc $outputdir/b1__b11-b001__b095-autosome-oxfordqc.fam
wc $outputdir/b1__b11-b001__b095-sexchrom-oxfordqc.fam

diff $outputdir/b1__b11-b001__b095-autosome-oxfordqc.fam $outputdir/b1__b11-b001__b095-sexchrom-oxfordqc.fam 

echo "Do we still have all the SNPs?"
wc $outputdir/b1__b11-b001__b095-autosome-oxfordqc.bim
wc $outputdir/b1__b11-b001__b095-sexchrom-oxfordqc.bim



# SPLIT BY CHROMOSOME
mkdir $outputdir/byChrom

for i in {1..22};
do \

echo $i

$plink --bfile $outputdir/b1__b11-b001__b095-autosome-oxfordqc --chr $i --keep-allele-order --make-bed --out $outputdir/byChrom/b1__b11-b001__b095-autosome-oxfordqc-chr$i

done;

for i in {23..26};
do \

echo $i

$plink --bfile $outputdir/b1__b11-b001__b095-sexchrom-oxfordqc --chr $i --keep-allele-order --make-bed --out $outputdir/byChrom/b1__b11-b001__b095-sexchrom-oxfordqc-chr$i

done;

# subset Y chromosome to just have males
mv $outputdir/byChrom/b1__b11-b001__b095-sexchrom-oxfordqc-chr24.bim $outputdir/byChrom/b1__b11-b001__b095-sexchrom-oxfordqc-chr24-all.bim
mv $outputdir/byChrom/b1__b11-b001__b095-sexchrom-oxfordqc-chr24.bed $outputdir/byChrom/b1__b11-b001__b095-sexchrom-oxfordqc-chr24-all.bed
mv $outputdir/byChrom/b1__b11-b001__b095-sexchrom-oxfordqc-chr24.fam $outputdir/byChrom/b1__b11-b001__b095-sexchrom-oxfordqc-chr24-all.fam

$plink --bfile $outputdir/byChrom/b1__b11-b001__b095-sexchrom-oxfordqc-chr24-all
$plink --bfile $outputdir/byChrom/b1__b11-b001__b095-sexchrom-oxfordqc-chr24-all --filter $outputdir/byChrom/b1__b11-b001__b095-sexchrom-oxfordqc-chr24-all.fam 1 -mfilter 3 --make-bed --out $outputdir/byChrom/b1__b11-b001__b095-sexchrom-oxfordqc-chr24


# WRITE SEX CHROMOSOMES WITH CHROMOSOME CODED 1, SO THAT BOLT WILL NOT FILTER THE SNPs!!

# XnonPAR
awk '$1 == 23 {print "1",$2,$3,$4,$5,$6} ' $outputdir/byChrom/b1__b11-b001__b095-sexchrom-oxfordqc-chr23.bim > $outputdir/byChrom/b1__b11-b001__b095-sexchrom-oxfordqc-chr23-recoded.bim

# Y
awk '$1 == 24 {print "1",$2,$3,$4,$5,$6} ' $outputdir/byChrom/b1__b11-b001__b095-sexchrom-oxfordqc-chr24.bim > $outputdir/byChrom/b1__b11-b001__b095-sexchrom-oxfordqc-chr24-recoded.bim

# PAR
awk '$1 == 25 {print "1",$2,$3,$4,$5,$6} ' $outputdir/byChrom/b1__b11-b001__b095-sexchrom-oxfordqc-chr25.bim > $outputdir/byChrom/b1__b11-b001__b095-sexchrom-oxfordqc-chr25-recoded.bim

# MT
awk '$1 == 26 {print "1",$2,$3,$4,$5,$6} ' $outputdir/byChrom/b1__b11-b001__b095-sexchrom-oxfordqc-chr26.bim > $outputdir/byChrom/b1__b11-b001__b095-sexchrom-oxfordqc-chr26-recoded.bim

