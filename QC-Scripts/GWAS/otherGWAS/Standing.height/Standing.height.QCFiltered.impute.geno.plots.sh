Rscript ../plot-GWAS-imputation-comparison.R /well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/otherGWAS/Standing.height/BOLTLMM.v16s/Standing.height-BOLT-LMM-v16s.out $SGE_TASK_ID plots -ymax 50 -qc -phenoFile /well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/otherGWAS/PhenotypesForBOLT-v3.txt -lreg -title -QCFiltered -bgenFile /well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/otherGWAS/Standing.height/BOLTLMM.v19s/Standing.height-BOLT-LMM-v19s-chr$SGE_TASK_ID.out -sample /well/ukbiobank/imputation/final/full/bgen/hrc+uk10k_sorted_8bit_rsids/imputed.sample/v2/chrPAR1.v2.sample -dontComputeLD -par 2