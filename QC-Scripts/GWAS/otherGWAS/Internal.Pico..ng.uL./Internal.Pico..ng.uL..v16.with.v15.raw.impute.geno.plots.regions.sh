Rscript ../plot-GWAS-imputation-comparison.R /well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/otherGWAS/Internal.Pico..ng.uL./BOLTLMM.v16/Internal.Pico..ng.uL.-BOLT-LMM-v16.out $SGE_TASK_ID plots -ymax 50 -phenoFile /well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/otherGWAS/PhenotypesForBOLT-v3.txt -lreg -title -raw -bgenFile /well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/otherGWAS/Internal.Pico..ng.uL./BOLTLMM.v15/Internal.Pico..ng.uL.-BOLT-LMM-v15-chr$SGE_TASK_ID.out -sample /well/ukbiobank/imputation/final/full/bgen/chr$SGE_TASK_ID.hrc+uk10k.I4.v1.1.sample -dontPlotManhattans -ldRData Internal.Pico..ng.uL.-BOLT-LMM-v16.out.chrgenome.maf0.001.miss0.05.pruned-QCFiltered-lreg-hitsAndLDcalc.RData -ldRDataImp Internal.Pico..ng.uL.-BOLT-LMM-v15-chr$SGE_TASK_ID.out.chr$SGE_TASK_ID.maf0.001.info0.3.pruned-QCFiltered-lreg-hitsAndLDcalc.RData
