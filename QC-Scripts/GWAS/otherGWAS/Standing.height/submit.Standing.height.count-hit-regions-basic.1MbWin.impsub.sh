#!/bin/bash
#$ -cwd
#$ -o Logs
#$ -j y
#$ -V

Rscript ../count-hit-regions-basic.R /well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/otherGWAS/Standing.height/BOLTLMM.v20b/Standing.height-BOLT-LMM-v20b-chr%%.out /well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/GWAS/otherGWAS/regions-hg19-window-1Mb.txt -title .1MbWin.IMPSUB -chr $SGE_TASK_ID -outdir /well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/otherGWAS/Standing.height/hitCounts
