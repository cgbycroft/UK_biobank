#!/bin/bash
#$ -cwd
#$ -o Logs
#$ -j y
#$ -V

Rscript ../compute-posteriors.R -gwas /well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/GWAS/otherGWAS/Standing.height/GIANT/GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq-with-positions.txt.gz -prior 0.2 -regions /well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/otherGWAS/Standing.height/hitCounts/regions-0.125cM-25KB-GIANT.chrgenome.txt -outdir /well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/otherGWAS/Standing.height/hitCounts -chr $SGE_TASK_ID -title .giantRegs
