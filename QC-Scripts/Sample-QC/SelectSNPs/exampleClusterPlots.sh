#################################
# Some cluster plots for talks etc
basedir=/well/ukbiobank/qcoutput.V2_QCed.sample-QC
cd $basedir

############
# Plate effect
# examples from:
# /well/ukbiobank/expt/V2_QCed.SNP-QC/data/SNP-QC_summary/V2_QCed.autosome.plate_effect.pval_lt_1e-12.UKBiLEVEAX_b1-b11.Batch_b001-b095.txt
Rscript /well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/R/scripts/V2_QCed.makeClusterPlots.R /well/ukbiobank/expt/V2_QCed.SNP-QC/data Affx-14796791,Affx-23780816,Affx-28637777,Affx-29504660 $basedir/QC-Scripts/Sample-QC/SelectSNPs/clusterPlots/plateEffect -b Batch_b001-Batch_b003 -colour Processed.Plate.Name


############
# Batch effect
# examples from:
# /well/ukbiobank/expt/V2_QCed.SNP-QC/data/SNP-QC_summary/V2_QCed.autosome.batch_effect.pval_lt_1e-12.UKBiLEVEAX_b1-b11.Batch_b001-b095.txt
Rscript /well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/R/scripts/V2_QCed.makeClusterPlots.R /well/ukbiobank/expt/V2_QCed.SNP-QC/data Affx-10361029,Affx-11485393,Affx-18681686,Affx-19597397 $basedir/QC-Scripts/Sample-QC/SelectSNPs/clusterPlots/batchEffect -b Batch_b001-Batch_b004


############
# Gender effect
# examples from:
# /well/ukbiobank/expt/V2_QCed.SNP-QC/data/SNP-QC_summary/V2_QCed.autosome.batch_effect.pval_lt_1e-12.UKBiLEVEAX_b1-b11.Batch_b001-b095.txt.pval_lt_1e-12.UKBiLEVEAX_b1-b11.Batch_b001-b095.txt
cut -f2 -d' ' /well/ukbiobank/expt/V2_QCed.SNP-QC/data/SNP-QC_summary/V2_QCed.autosome.malefemale.pval_lt_1e-12.UKBiLEVEAX_b1-b11.Batch_b001-b095.txt | sort -u > $basedir/QC-Scripts/Sample-QC/SelectSNPs/clusterPlots/genderEffect/gender.effect.autosome.snps.txt
# just plot all 8 autosomal snps
Rscript /well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/R/scripts/V2_QCed.makeClusterPlots.R /well/ukbiobank/expt/V2_QCed.SNP-QC/data $basedir/QC-Scripts/Sample-QC/SelectSNPs/clusterPlots/genderEffect/gender.effect.autosome.snps.txt $basedir/QC-Scripts/Sample-QC/SelectSNPs/clusterPlots/genderEffect -b Batch_b001-Batch_b003 -sex


############
# Array effect
# examples from:
# /well/ukbiobank/expt/V2_QCed.SNP-QC/data/SNP-QC_summary/V2_QCed.autosome.array_test.pval_lt_1e-12.UKBiLEVEAX_b1-b11.Batch_b001-b095.txt.pval_lt_1e-12.UKBiLEVEAX_b1-b11.Batch_b001-b095.txt
Rscript /well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/R/scripts/V2_QCed.makeClusterPlots.R /well/ukbiobank/expt/V2_QCed.SNP-QC/data Affx-10096426,Affx-10115878,Affx-10382633,Affx-10412935 $basedir/QC-Scripts/Sample-QC/SelectSNPs/clusterPlots/arrayEffect -b UKBiLEVEAX_b1-UKBiLEVEAX_b9
Rscript /well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/R/scripts/V2_QCed.makeClusterPlots.R /well/ukbiobank/expt/V2_QCed.SNP-QC/data Affx-10096426,Affx-10115878,Affx-10382633,Affx-10412935 $basedir/QC-Scripts/Sample-QC/SelectSNPs/clusterPlots/arrayEffect -b Batch_b001-Batch_b009


############
# HWE
# examples from:
# /well/ukbiobank/expt/V2_QCed.SNP-QC/data/SNP-QC_summary/V2_QCed.autosome.hwe.pval_lt_1e-12.UKBiLEVEAX_b1-b11.Batch_b001-b095.txt
Rscript /well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/R/scripts/V2_QCed.makeClusterPlots.R /well/ukbiobank/expt/V2_QCed.SNP-QC/data Affx-12351817,Affx-12376960,Affx-15915035 $basedir/QC-Scripts/Sample-QC/SelectSNPs/clusterPlots/hweTest -b UKBiLEVEAX_b1-UKBiLEVEAX_b9
Rscript /well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/R/scripts/V2_QCed.makeClusterPlots.R /well/ukbiobank/expt/V2_QCed.SNP-QC/data Affx-10181029,Affx-11934703,Affx-14343512 $basedir/QC-Scripts/Sample-QC/SelectSNPs/clusterPlots/hweTest -b Batch_b001-Batch_b009




############
# Array effect test --> EXTRA FROM COLIN Feb 2017. Failed in Male or Females, but not in combined.
# examples from:
# /well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/Sample-QC/SelectSNPs/ExtraSNPsArrayFail-inSampleQC.txt
Rscript /well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/R/scripts/V2_QCed.makeClusterPlots.R /well/ukbiobank/expt/V2_QCed.SNP-QC/data Affx-34704662,Affx-34574387,Affx-34849471 $basedir/QC-Scripts/Sample-QC/SelectSNPs/clusterPlots/ArrayTestSexInSampleQC -b UKBiLEVEAX_b1,UKBiLEVEAX_b5,UKBiLEVEAX_b11,Batch_b001,Batch_b059,Batch_b092
# CONCLUSION ==> Most likely genotype calling errors, rather than problems with intensities ==> they won't badly effect the L2R statistics.
