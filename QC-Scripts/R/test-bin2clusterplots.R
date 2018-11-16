

## I plan to create a package out of the R scripts...
source("/well/ukbiobank/qcoutput/QC-Scripts/R/scripts/bin2clusterplots.R")
source("/well/ukbiobank/qcoutput/QC-Scripts/R/scripts/batch2sampleinfo.R")


Batch = "Batch_b080"
BatchNo = "b080"

## A directory which contains :
## AxiomGT1.autosome.calls.bin, AxiomGT1.autosome.intensities.bin, AxiomGT1.autosome.snp-posteriors.bin
## AxiomGT1.sexchrom.calls.bin, AxiomGT1.sexchrom.intensities.bin, AxiomGT1.sexchrom.snp-posteriors.bin
datadir = paste0("/well/ukbiobank/qcoutput/",Batch,"/GT1")

## The location and name of the augmented sample csv table
## Will be read in to get information about gender, ethnicity, plate, etc.
phenodir = paste0("/well/ukbiobank/qcoutput/",Batch)
phenofile = paste0("UKB_WCSGAX_",BatchNo,"_Sample_Table_Pheno.csv")



## The basic cluster plot


## The index of the probeset to be plotted, i.e.,
## the row number in the binary tables of calls and intensities
pid = 1 

## Is the probeset on chromosomes 1..22 or not? It is necessary to specify this information
## becase the intensities are stored into two different files, for each batch.
is.autosomal = TRUE

## It might be useful to check that the row number (pid) corresponds to the probeset you
## have in mind (but pname can left NULL as well)
## You can look up the list of probesets (in the order they are in the binary files) here:
## /well/ukbiobank/qcoutput/QC-Scripts/Axiom_UKBL_Chip (for UK BiLEVE batches)
## /well/ukbiobank/qcoutput/QC-Scripts/Axiom_UKBB_Chip (for UK Biobank batches)

## If the probeset name in the binary files does not much the optional pname argument, then NULL is returned
pname = "XXXXXXXXXXXXXXX"
data = unpack.probeset.data(datadir,pid=pid,is.autosomal=is.autosomal,pname=pname)

## Otherwise, intensity data is returned
pname = "AFFX-KIT-000001"
data = unpack.probeset.data(datadir,pid=pid,is.autosomal=is.autosomal,pname=pname)


png(file="clusterplot-example-autosomal.png",
    width=7,height=7,units="in",res=150)
cluster.plot(data,main.text="An informative title")
dev.off( )


## Works much the same way with the sex-linked probesets

data = unpack.probeset.data(datadir,pid=1,is.autosomal=FALSE)
png(file="clusterplot-example-sexchrom.png",
    width=7,height=7,units="in",res=150)
cluster.plot(data,main.text="An informative title")
dev.off( )



## By default the points in a cluster plots are colored by genotype
## but it is easy to color them according to some other categorization
## First read in the sample information
CsvFile = paste(phenodir,"/",phenofile,sep="")
BatchInfo = get.batch.info(Batch,CsvFile)


## Assign blue to males and red to females
Colors = c("blue","red")
names(Colors) = c("M","F")


## Cluster plots of sex-linked SNPs make more sense if
## males and females are colored differently
png(file="clusterplot-example-sexchrom-color-by-gender.png",
    width=7,height=7,units="in",res=150)
cluster.plot(data,main.text="An informative title",
	     col=Colors[BatchInfo$Inferred.Gender])
dev.off( )
## And here is how to color by self-reported ethinicity
png(file="clusterplot-example-sexchrom-color-by-pop.png",
    width=7,height=7,units="in",res=150)
cluster.plot(data,main.text="An informative title",
	     col=BatchInfo$Colors,
	     pch=BatchInfo$Chars)
dev.off( )
