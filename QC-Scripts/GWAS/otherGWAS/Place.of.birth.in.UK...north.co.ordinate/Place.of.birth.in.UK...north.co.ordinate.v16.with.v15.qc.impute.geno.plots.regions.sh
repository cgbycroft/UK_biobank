Rscript ../plot-GWAS-imputation-comparison.R /well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/otherGWAS/Place.of.birth.in.UK...north.co.ordinate/BOLTLMM.v16/Place.of.birth.in.UK...north.co.ordinate-BOLT-LMM-v16.out $SGE_TASK_ID plots -ymax 50 -phenoFile /well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/otherGWAS/PhenotypesForBOLT-v3.txt -qc -lreg -title -QCFiltered -bgenFile /well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/otherGWAS/Place.of.birth.in.UK...north.co.ordinate/BOLTLMM.v15/Place.of.birth.in.UK...north.co.ordinate-BOLT-LMM-v15-chr$SGE_TASK_ID.out -dontPlotManhattans -ldRData Place.of.birth.in.UK...north.co.ordinate-BOLT-LMM-v16.out.chrgenome.maf0.001.miss0.05.pruned-QCFiltered-lreg-hitsAndLDcalc.RData -ldRDataImp Place.of.birth.in.UK...north.co.ordinate-BOLT-LMM-v15-chr$SGE_TASK_ID.out.chr$SGE_TASK_ID.maf0.001.info0.3.pruned-QCFiltered-lreg-hitsAndLDcalc.RData
