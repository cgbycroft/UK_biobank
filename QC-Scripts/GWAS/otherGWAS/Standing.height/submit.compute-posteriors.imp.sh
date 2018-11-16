#!/bin/bash
#$ -cwd
#$ -o Logs
#$ -j y
#$ -V

Rscript ../compute-posteriors.R -gwas /well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/otherGWAS/Standing.height/BOLTLMM.v20/Standing.height-BOLT-LMM-v20-chr%%.out -prior 0.2 -regions /well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/otherGWAS/Standing.height/hitCounts/regions-0.125cM-25KB-GIANT.chrgenome.txt -outdir /well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/otherGWAS/Standing.height/hitCounts -chr $SGE_TASK_ID -title .giantRegs
