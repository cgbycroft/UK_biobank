

## I plan to create a package out of the R scripts...
source("/well/ukbiobank/qcoutput/QC-Scripts/R/scripts/bin2clusterplots.R")
source("/well/ukbiobank/qcoutput/QC-Scripts/R/scripts/batch2sampleinfo.R")


Batch = "Batch_b081"
BatchNo = "b081"


## A directory which contains :
## AxiomGT1.autosome.calls.bin, AxiomGT1.autosome.intensities.bin, AxiomGT1.autosome.snp-posteriors.bin
## AxiomGT1.sexchrom.calls.bin, AxiomGT1.sexchrom.intensities.bin, AxiomGT1.sexchrom.snp-posteriors.bin
datadir = paste0("/well/ukbiobank/qcoutput/",Batch,"/GT1")

## The location and name of the augmented sample csv table
## Will be read in to get information about gender, ethnicity, plate, etc.
phenodir = paste0("/well/ukbiobank/qcoutput/",Batch)
phenofile = paste0("UKB_WCSGAX_",BatchNo,"_Sample_Table_Pheno.csv")



## Is the probeset on chromosomes 1..22 or not? It is necessary to specify this information
## becase the intensities are stored into two different files, for each batch.
is.autosomal = TRUE

## It might be useful to check that the row number (pid) corresponds to the probeset you
## have in mind (but pname can left NULL as well)
## You can look up the list of probesets (in the order they are in the binary files) here:
## /well/ukbiobank/qcoutput/QC-Scripts/Axiom_UKBL_Chip (for UK BiLEVE batches)
## /well/ukbiobank/qcoutput/QC-Scripts/Axiom_UKBB_Chip (for UK Biobank batches)

## Cluster plot with colors
CsvFile = paste(phenodir,"/",phenofile,sep="")
BatchInfo = get.batch.info(Batch,CsvFile)


## Color males in blue and females in red
Colors = c("blue","red")
names(Colors) = c("M","F")


png(file=paste0(Batch,"-SNPs-fail-Males2Females-autosome%02d.png"),
    width=7,height=7,units="in",res=150)

## These probesets turn out to have very small p-values for the
## SNP QC test which compares the genotype frequencies of males
## to those of females
for (pid in c(171019,
	      423165,
	      435933,
	      310168,
	      357151,
	      317201,
	      673704)) {

  data = unpack.probeset.data(datadir,pid=pid,is.autosomal=is.autosomal)

  ## Cluster plots of sex-linked SNPs make more sense if
  ## males and females are colored differently
  cluster.plot(data,main.text="An informative title",
	       col=Colors[BatchInfo$Inferred.Gender]) 
}

dev.off( )
