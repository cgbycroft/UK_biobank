plink=/well/ukbiobank/qcoutput/Src/plink-v1.90b3k-x86_64/plink

############ %%%
inputdir=$basedir/data/Combined/byChrom
############

# Oxford qc'd files
for i in `seq 1 22`;
do \
    
echo $i    
$plink --bfile $inputdir/b1__b11-b001__b095-autosome-oxfordqc-chr$i --keep-allele-order --freq --hardy --het gz --missing gz --out $inputdir/b1__b11-b001__b095-autosome-oxfordqc-chr$i

done;
