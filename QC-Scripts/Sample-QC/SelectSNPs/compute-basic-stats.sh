################################
# For UKBiobank sample QC pipeline step 1
################################
# Clare Bycroft, May 2016
#
# This script computes basic statistics using plink, on the Combined plink files, after running merge-batches.sh This is before any filtering is done with regards to relatedness etc.
################################

plink=/well/ukbiobank/qcoutput/Src/plink-v1.90b3k-x86_64/plink

############
inputdir=$basedir/data/ByBatch # this is the directory with plink files made by create-plink-subsets.sh
outputdir=$basedir/data/Combined
############

# AUTOSOMES

# freq; hwe; het; missing
$plink --bfile $outputdir/b1__b11-b001__b095-autosome-sampleqc --keep-allele-order --freq --hardy --het gz --missing gz --out $outputdir/b1__b11-b001__b095-autosome-sampleqc


# SEX CHROMOSOMES

# freq; hwe; het; missing
$plink --bfile $outputdir/b1__b11-b001__b095-sexchrom-sampleqc --keep-allele-order --freq --hardy --het gz --missing gz --out $outputdir/b1__b11-b001__b095-sexchrom-sampleqc


# Different sex chroms separately
for i in {23..26};
do \
  
    $plink --bfile $outputdir/b1__b11-b001__b095-sexchrom-sampleqc --keep-allele-order --chr $i --freq --hardy --missing gz --make-bed --out $outputdir/b1__b11-b001__b095-sexchrom-sampleqc-$i   

    # for heterozygosity, need to trick plink
    awk '{ print "1",$2,$3,$4,$5,$6'} $outputdir/b1__b11-b001__b095-sexchrom-sampleqc-$i.bim > tmp-$i.bim    

    $plink --bim tmp-$i.bim --bed $outputdir/b1__b11-b001__b095-sexchrom-sampleqc-$i.bed --fam $outputdir/b1__b11-b001__b095-sexchrom-sampleqc-$i.fam --keep-allele-order --het gz --out $outputdir/b1__b11-b001__b095-sexchrom-sampleqc-$i
    
done;
