#!/bin/bash
#$ -cwd
#$ -o Logs
#$ -j y
#$ -V

Rscript ../count-hit-regions-basic.R /well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/GWAS/otherGWAS/Standing.height/GIANT/GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq-with-positions.txt.gz /well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/GWAS/otherGWAS/regions-hg19-window-1Mb.txt -title .1MbWin.GIANT -chr $SGE_TASK_ID -outdir /well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/otherGWAS/Standing.height/hitCounts
