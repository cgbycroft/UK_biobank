#!/bin/bash
#$ -cwd
#$ -o Logs
#$ -j y
#$ -V

Rscript ../compute-bayes-factors.R -gwas /well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/otherGWAS/Standing.height/BOLTLMM.v20/Standing.height-BOLT-LMM-v20-chr%%.out -chr $SGE_TASK_ID -prior 0.2
