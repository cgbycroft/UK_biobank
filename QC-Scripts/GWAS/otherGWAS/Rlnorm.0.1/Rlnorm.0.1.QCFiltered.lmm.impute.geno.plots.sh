Rscript ../plot-GWAS-imputation-comparison.R /well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/otherGWAS/Rlnorm.0.1/BOLTLMM.v18/Rlnorm.0.1-BOLT-LMM-v18.out $SGE_TASK_ID plots -ymax 50 -qc -phenoFile /well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/otherGWAS/PhenotypesForBOLT-v3.txt -title -QCFiltered -bgenFile /well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/otherGWAS/Rlnorm.0.1/BOLTLMM.v17/Rlnorm.0.1-BOLT-LMM-v17-chr$SGE_TASK_ID.out -sample /well/ukbiobank/imputation/final/full/bgen/hrc+uk10k_sorted_v1.1/chr$SGE_TASK_ID.hrc+uk10k.I4.v1.1.sample -dontComputeLD
