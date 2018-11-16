#################################
# This file contains descriptions of each of the fields in the sample table (genotypes) that will be provided to researchers.
# This file is read in by create-sample-table.R to create a new file that corresponds to a particular version of the table.
#################################

Best.Array <- c("string","Affymetrix string which identifes a .CEL file, and also corresponds uniquely to a DNA sample in the released genotype calls files.")
Sample.Name <- c("string","ID created by UKBiobank for each DNA sample extracted. This corresponds to a single DNA sample in the released genotype calls files. This id is only unique within each of the subsets typed on the two arrays (UKBiobank and UKBiLEVE)")
genotyping.array <- c("string","Indicates which array the DNA sample was typed on. UKBL=UKBiLEVE and UKBB=UKBiobank")
Batch <- c("string","Indicates in which batch Affymetrix called the genotypes for this sample. See Affymetrix documentation for details.")
Plate.Name <- c("string","Indicates on which plate Affymetrix processed this sample. See Affymetrix documentation for details.")
Well <- c("string","Indicates in which well Affymetrix processed sample. See Affymetrix documentation for details.")
Cluster.CR <- c("numeric","Affymetrix quality control metric. Corresponds to sample-specific call rate for each individual, computed using probesets internal to Affymetrix. See Affymetrix SNP-polisher documentation for details.")
dQC <-c("numeric","Affymetrix quality control metric. Measures the resolution of the distributions of intensity 'contrast' values, based on intensities of probe sequences for non-polymorphic genome locations. See Affymetrix SNP-polisher documentation for details. Note that values are rounded to 2 decimal places for individuals type on the UKBibank array, and 5 decimal places for individuals typed on the UKBiLEVE array.")
Internal.Pico..ng.uL. <- c("numeric","Affymetrix quality control metric. Measure of DNA concentration in sample. See Affymetrix documentation for details.")
Submitted.Gender <- c("string","Gender submitted by participant. Should match UKBiobank phenotype field 'Sex'.")
Inferred.Gender <- c("string","Sex as determined by Affymetrix, and used in calling genotypes on the Y chromosome and the sex-specific region of the X chromosome. See Affymetrix documentation for details.")
X.intensity <- c("numeric","Affymetrix metric used for determining sex. Measures average probe intensity on a set of non-polymorphic probes on the X chromosome. Note this set of probes are not in the genotype calls.")
Y.intensity <- c("numeric","Affymetrix metric used for determining sex. Measures average probe intensity on a set of non-polymorphic probes on the Y chromosome. Note this set of probes are not in the genotype calls.")
Submitted.Plate.Name <- c("string","Indicates on which plate the DNA sample was delivered to Affymetrix. See UKBiobank DNA processing documentation for details.")
Submitted.Well <- c("string","Indicates in which well the DNA sample was delivered to Affymetrix. See UKBiobank DNA processing documentation for details.")
sample.qc.missing.rate <-  c("numeric","Missing rate of each sample based on a set of high-quality SNPs. This metric was used to identify outliers based on heterozygosity and missing rates (het.missing.outliers). See genotype QC documentation for details.")
heterozygosity <- c("numeric","Heterozygosity across a set of high-quality SNPs. See genotype QC documentation for details.")
heterozygosity.pc.corrected <- c("numeric","Heterozygosity after adjusting for ancestry using principal components. This metric was used to identify outliers based on heterozygosity and missing rates (het.missing.outliers). See genotype QC documentation for details.")
het.missing.outliers <- c("indicator","Indicates samples identified as outliers in heterozygosity and missing rates, which indicates poor-quality genotypes for these samples.")
putative.sex.chromosome.aneuploidy <- c("indicator","Indicates samples identified as putatively carrying sex chromosome configurations that are not either XX or XY. These were intified by looking at average log2Ratios for Y and X chromosomes. See genotype QC documentation for details.")
in.kinship.table <- c("indicator","Indicates samples which have at least one relative among the set of genotyped individuals. These are exactly the set of samples that appear in the kinship table. See genotype QC documentation for details.")
excluded.from.kinship.inference <- c("indicator","Indicates samples which were excluded from the kinship inference procedure. See genotype QC documentation for details.")
excess.relatives <- c("indicator","Indicates samples which have more than 10 putative third-degree relatives in the kinship table.")
in.white.British.ancestry.subset <-  c("indicator","Indicates samples who self-reported 'White British' and have very similar genetic ancestry based on a principal components analysis of the genotypes. See genotype QC documentation for details.")
used.in.pca.calculation <- c("indicator","Indicates samples which were used to compute principal components. All samples with genotype data were subsequently projected on to each of the 40 computed components. See genotype QC documentation for details.")
PC1 <- c("numeric","Score on principal component 1. See genotype QC documentation for details.")
PC2 <- c("numeric","Score on principal component 2. See genotype QC documentation for details.")


# Kinship table fields
ID1 <- c("string","Best.Array string for individual 1 in related pair.")
ID2 <- c("string","Best.Array string for individual 2 in related pair.")
HetHet <- c("numeric","Fraction of SNPs for which the pair both have a heterozygous genotype (output from KING software).")
IBS0 <- c("numeric","Fraction of SNPs for which the pair shares zero alleles (output from KING software).")
Kinship <- c("numeric","Estimate of the kinship coefficient for this pair based on the set of SNPs used in the kinship inference (output from KING software). The set of SNPs is indicated by the field: used.in.kinship.inference")
