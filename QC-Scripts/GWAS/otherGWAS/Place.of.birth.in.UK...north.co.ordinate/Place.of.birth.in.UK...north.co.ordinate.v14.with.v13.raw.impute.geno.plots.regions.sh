Rscript ../plot-GWAS-imputation-comparison.R /well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/otherGWAS/Place.of.birth.in.UK...north.co.ordinate/BOLTLMM.v14/Place.of.birth.in.UK...north.co.ordinate-BOLT-LMM-v14.out $SGE_TASK_ID plots -ymax 50 -phenoFile /well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/otherGWAS/PhenotypesForBOLT-v3.txt -lreg -title -raw -bgenFile /well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/otherGWAS/Place.of.birth.in.UK...north.co.ordinate/BOLTLMM.v13/Place.of.birth.in.UK...north.co.ordinate-BOLT-LMM-v13-chr$SGE_TASK_ID.out -dontPlotManhattans -ldRData Place.of.birth.in.UK...north.co.ordinate-BOLT-LMM-v14.out.chrgenome.maf0.001.miss0.05.pruned-QCFiltered-lreg-hitsAndLDcalc.RData