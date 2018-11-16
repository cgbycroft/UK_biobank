##################
# Script to extract allele frequencies for Europeans from large vcf files.
##################

# cd /well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/Exac
# wget ftp://ftp.broadinstitute.org/pub/ExAC_release/release1/ExAC.r1.sites.vep.vcf.gz
# wget ftp://ftp.broadinstitute.org/pub/ExAC_release/release1/ExAC.r1.sites.vep.vcf.gz.tbi

vcftools=/apps/well/vcftools/0.1.14-gcc4.7.2/bin/vcftools
vcfsubset=/apps/well/vcftools/0.1.14-gcc4.7.2/bin/vcf-subset
bcftools=/apps/well/bcftools/1.3/bin/bcftools
plink=/well/ukbiobank/qcoutput/Src/plink-v1.90b3k-x86_64/plink


basedir=/well/ukbiobank/qcoutput.V2_QCed.sample-QC
outdir=$basedir/data/Exac
exac=/well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/Exac/ExAC.r1.sites.vep.vcf.gz
ukbbgeno=/well/ukbiobank/expt/V2_QCed.export/data/calls_export/v3/V2_QCed.export.UKBiLEVEAX_b1-b11.Batch_b001_b095


echo "using data from..."
echo $ukbbgeno
echo $exac


# mkdir $outdir

cd $basedir/QC-Scripts/ExacComparison

	
######## Get list of UKBB snps
awk 'OFS="\t" {print $1,$4}' $ukbbgeno*.bim > ukbbSNPs.txt


######## Extract ExAC snps and allele count columns using bcftools.
# Takes about half an hour.?

# AC_NFE = Alternative allele count non-Finnish Europeans.
# AN_NFE = Non-missing alleles count non-Finnish Europeans.
# zcat $exac | grep '##INFO' 

#CHROM	POS	ID	REF	ALT	QUAL	FILTER

time $bcftools view --regions-file ukbbSNPs.txt -O z --output-file $outdir/ExAC.r1.sites.vep.UKBB-SNPs.vcf.gz $exac

# exctract suitable columns (fast once subsetting has been done)
time $bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AC_NFE\t%INFO/AN_NFE\n' --output-file $outdir/ExAC.r1.sites.vep.UKBB-SNPs.aCounts_NFE.txt $outdir/ExAC.r1.sites.vep.UKBB-SNPs.vcf.gz

