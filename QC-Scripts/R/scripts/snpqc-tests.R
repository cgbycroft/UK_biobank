

source("/well/ukbiobank/qcoutput.V2_QCed.sample-QC-testing/QC-Scripts/R/scripts/batch2sampleinfo.R")
library(dplyr)
library(UKBiobankSubroutines)

pidSize = 16

UKBBChipDir = "/well/ukbiobank/expt/V2_QCed.SNP-QC/data/V2_QCed.Axiom_UKBB_Chip"
UKBLChipDir = "/well/ukbiobank/expt/V2_QCed.SNP-QC/data/V2_QCed.Axiom_UKBL_Chip"


## Plink uses the following codes to specify non-autosome chromosome types:
## X  -> 23 (X chromosome)
## Y  -> 24 (Y chromosome)
## XY -> 25 (Pseudo-autosomal region)
## MT -> 26 (Mitochondrial)
chrom.as.numeric <- function(ps2snp) {
    Chrom = as.character(ps2snp$Chromosome)
    Chrom = factor(Chrom,levels=c(1:22,"X","Y","XY","MT","---"),labels=1:27)
    Chrom = as.numeric(as.character(Chrom))
    IsPAR = (ps2snp$IsInPAR1+ps2snp$IsInPAR2)>0
    Chrom[IsPAR] = 25
    return(Chrom)
}
ukbileve.ps2snp <- function(ProbeSet) {
    ps2snp = read.table(paste(UKBLChipDir,"/",ProbeSet,".txt",sep=""),header=TRUE,stringsAsFactors=FALSE)
    ps2snp$Chromosome = chrom.as.numeric(ps2snp)
    return(ps2snp)
}
ukbiobank.ps2snp <- function(ProbeSet) {
    ps2snp = read.table(paste(UKBBChipDir,"/",ProbeSet,".txt",sep=""),header=TRUE,stringsAsFactors=FALSE)
    ps2snp$Chromosome = chrom.as.numeric(ps2snp)
    return(ps2snp)
}
ukb.axiom.AffySNPID <- function(ProbeSet) {
  UKBL.ps2snp = read.table(paste(UKBLChipDir,"/",ProbeSet,".txt",sep=""),header=TRUE,stringsAsFactors=FALSE)
  UKBB.ps2snp = read.table(paste(UKBBChipDir,"/",ProbeSet,".txt",sep=""),header=TRUE,stringsAsFactors=FALSE)
  UKBL.ps2snp = tbl_df(UKBL.ps2snp)
  UKBB.ps2snp = tbl_df(UKBB.ps2snp)
  UKBL.ps2snp = select(UKBL.ps2snp, AffySNPID, Chromosome, Position, IsInPAR1, IsInPAR2)
  UKBB.ps2snp = select(UKBB.ps2snp, AffySNPID, Chromosome, Position, IsInPAR1, IsInPAR2)
  ps2snp = merge(UKBL.ps2snp, UKBB.ps2snp, all = TRUE)
  ps2snp = ps2snp[!duplicated(ps2snp, fromLast = TRUE),]
  ps2snp$Chromosome = chrom.as.numeric(ps2snp)
  return(ps2snp)
}
ukb.axiom.ProbeSetID <- function(ProbeSet) {
  UKBL.ps2snp = read.table(paste(UKBLChipDir,"/",ProbeSet,".txt",sep=""),header=TRUE,stringsAsFactors=FALSE)
  UKBB.ps2snp = read.table(paste(UKBBChipDir,"/",ProbeSet,".txt",sep=""),header=TRUE,stringsAsFactors=FALSE)
  UKBL.ps2snp = tbl_df(UKBL.ps2snp)
  UKBB.ps2snp = tbl_df(UKBB.ps2snp)
  UKBL.ps2snp = select(UKBL.ps2snp, ProbeSetID, Chromosome, Position, IsInPAR1, IsInPAR2)
  UKBB.ps2snp = select(UKBB.ps2snp, ProbeSetID, Chromosome, Position, IsInPAR1, IsInPAR2)
  ps2snp = merge(UKBL.ps2snp, UKBB.ps2snp, all = TRUE)
  ps2snp = ps2snp[!duplicated(ps2snp, fromLast = TRUE),]
  ps2snp$Chromosome = chrom.as.numeric(ps2snp)
  return(ps2snp)
}

SummaryCountsFromAxiomCallsFile <- function(ProbeSet,ps2snp,
					    BinDatadir,BatchInfo,
					    PsPerformance,CountsFile,
					    exclude.indivs=character()) {
  writeLines("SummaryCountsFromAxiomCallsFile :")
  if (sub.is.ukbiobank(Batch)) {
    ProbeSetID = ukbiobank.ps2snp(ProbeSet)$ProbeSetID
  } else if (sub.is.ukbileve(Batch)) {
    ProbeSetID = ukbileve.ps2snp(ProbeSet)$ProbeSetID
  } else {
    stop(paste("Batch ",Batch," not found"))
  }

  IDsToExclude = match(exclude.indivs,BatchInfo$PIID)
  
  CallsConnection = file(paste(BinDatadir,"/AxiomGT1.",ProbeSet,".calls.bin",sep=""),"rb")
  values = readBin(CallsConnection,numeric(),n=3)
  nProbes = as.integer(values[1])
  nSamples = as.integer(values[2])
  recSize = as.integer(values[3])
  
  writeLines(c("ProbeSet = ",ProbeSet))
  writeLines(c("nSamples = ",nSamples))
  writeLines(c("nProbesets = ",nProbes))
  fromProbe = 1
  toProbe = nProbes
  
  Gender = BatchInfo$Inferred.Gender
  if (length(Gender)!=nSamples) { stop("length(Gender)!=nSample") }
  if (!is.factor(Gender)) { Gender = as.factor(Gender) }
  
  ## Read the data in chunks of size = chunkSize
  probesToRead = toProbe - fromProbe + 1
  chunkSize = 10000
  firstChunk = TRUE
  
  U = "HG00097"
  V = "HG00264"
  uReferences = get.unique.references(BatchInfo,U)
  vReferences = get.unique.references(BatchInfo,V)
  allReferences = get.batch.references(BatchInfo)
  IncludeIndivs = 1:nSamples
  IncludeIndivs = setdiff(IncludeIndivs,allReferences)            ## always exclude the reference samples
  IncludeIndivs.AF = intersect(IncludeIndivs,which(Gender=='F'))  ## all_female
  IncludeIndivs.AM = intersect(IncludeIndivs,which(Gender=='M'))  ## all_male
  IncludeIndivs = setdiff(IncludeIndivs,IDsToExclude)             ## exclude outliers
  IncludeIndivs.WF = intersect(IncludeIndivs,which(Gender=='F'))  ## british_white_female
  IncludeIndivs.WM = intersect(IncludeIndivs,which(Gender=='M'))  ## british_white_male
  
  nFemales = length(IncludeIndivs.WF)
  nMales = length(IncludeIndivs.WM)
  writeLines(c("Number of females = ",nFemales))
  writeLines(c("Number of males = ",nMales))
  
  if (fromProbe>1) {
    seek(CallsConnection,(fromProbe-1)*(1*nSamples+pidSize),origin="current")
  } else {
    Clusters = c("AA","AB","BB","NC")
    ColumnsNames = c("ProbeSetID",
      paste(Clusters,"all_female",sep="_"),paste(Clusters,"all_male",sep="_"),
      paste(Clusters,"ceu_female",sep="_"),paste(Clusters,"ceu_male",sep="_"),
      paste(Clusters,   "HG00097",sep="_"),paste(Clusters, "HG00264",sep="_"))
  }
  
  Counts = data.frame( )
  writeLines("Processing...")
  while (probesToRead>0) {
    
    chunkSize = min(chunkSize,probesToRead)
    
    Bytes = readBin(CallsConnection,what="raw",n=chunkSize*(1*nSamples+pidSize))
    CharBytes = rep(1:pidSize,times=chunkSize) + rep(0:(chunkSize-1),each=pidSize)*(1*nSamples+pidSize)
    ProbeSets = readBin(Bytes[CharBytes],what="character",n=chunkSize)
    Calls = readBin(Bytes[-CharBytes],what="int",size=1,n=chunkSize*nSamples)
    Calls = matrix(Calls,ncol=nSamples,byrow=TRUE)
    
    chunkProbeSetID = ProbeSetID[1:chunkSize]
    ProbeSetID = ProbeSetID[-(1:chunkSize)]
    ProbeSets = sub(" *","",ProbeSets)
    
    if(sum(ProbeSets!=chunkProbeSetID)) {
      stop("ProbeSets!=chunkProbeSetID")
    }
    if (firstChunk) {
      writeLines("/* Small chunk of the calls matrix */")
      print(Calls[1:4,1:8])
    }
    
    Counts.AF = genotype.counts(Calls[,IncludeIndivs.AF])
    Counts.AM = genotype.counts(Calls[,IncludeIndivs.AM])
    Counts.WF = genotype.counts(Calls[,IncludeIndivs.WF])
    Counts.WM = genotype.counts(Calls[,IncludeIndivs.WM])
    uCounts = genotype.counts(Calls[,uReferences])
    vCounts = genotype.counts(Calls[,vReferences])
    
    Chunk = data.frame(ProbeSets,Counts.AF,Counts.AM,Counts.WF,Counts.WM,uCounts,vCounts,stringsAsFactors = FALSE)
    Counts = rbind(Counts,Chunk)
    
    firstChunk = FALSE
    probesToRead = probesToRead - chunkSize
    print(toProbe - probesToRead,quote=FALSE)
  }
  
  colnames(Counts) = ColumnsNames
  close(CallsConnection)

  ApplyAffymetrixQCtoAffySNPs(ps2snp,PsPerformance,Counts,CountsFile)
}
filteredFromPs.performance <- function(PsPerformance) {

  filtered = read.table(PsPerformance,header=TRUE,stringsAsFactors=FALSE)
  filtered = tbl_df(filtered)

  if (match("affy_snp_id",names(filtered),nomatch=0)) {
    filtered = rename(filtered,snpid=affy_snp_id)
  }
  if (!match("probeset_id",names(filtered),nomatch=0) ||
      !match("snpid",names(filtered),nomatch=0) ||
      !match("ConversionType",names(filtered),nomatch=0) ||
      !match("BestProbeset",names(filtered),nomatch=0) ) {
    print(head(filtered))
    stop("head(filtered) : no probeset_id, or snpid, or ConversionType, or BestProbeset")
  }
  filtered = select(filtered, probeset_id, snpid, ConversionType, BestProbeset)

  ## Filter based on whether the probeset is the best probeset
  filtered = filter(filtered, BestProbeset == 1)
  ## Filter based on whether the probeset is converted or not
  filtered = filter(filtered, ConversionType == "PolyHighResolution" | ConversionType == "NoMinorHom" |
		      ConversionType == "MonoHighResolution" | ConversionType == "Hemizygous")

  return(filtered)
}
ApplyAffymetrixQCtoAffySNPs <- function(ps2snp,PsPerformance,Data,OutputFile) {
  filtered = filteredFromPs.performance(PsPerformance)

  ## Join the data and the SNP information (at filtered probesets only)
  filteredData = inner_join(filtered, Data, by = c("probeset_id" = "ProbeSetID"))
  ## The probeset_id, ConversionType and BestProbeset are no longer needed
  filteredData = select(filteredData,-probeset_id,-ConversionType,-BestProbeset)

  ## Align the SNPs with actual data with all the SNPs
  ## The SNPs with missing data get assigned NAs
  AffySNPID = data.frame(AffySNPID = ps2snp$AffySNPID, stringsAsFactors = FALSE)
  filteredData = left_join(AffySNPID,filteredData, by = c("AffySNPID" = "snpid"))

  Connection = file(OutputFile,"wt")
  filteredData = format(filteredData,scientific=TRUE,digits=9)
  write.table(filteredData,Connection,quote=FALSE,row.names=FALSE,col.names=TRUE)
  close(Connection)
}
AggregateCountsFromCountsFiles <- function(ProbeSet,ps2snp,CountsFiles,TotalFile) {
  writeLines("AggregateCountsFromCountsFiles : ")
  ProbeSetID = "AffySNPID"
  
  Clusters = c("AA","AB","BB","NC")
  RowNames = ps2snp[[ProbeSetID]]
  ColNames = c(ProbeSetID,
    paste(Clusters,"all_female",sep="_"),paste(Clusters,"all_male",sep="_"),
    paste(Clusters,"ceu_female",sep="_"),paste(Clusters,"ceu_male",sep="_"),
    paste(Clusters,   "HG00097",sep="_"),paste(Clusters, "HG00264",sep="_"))
  nRows = length(RowNames)
  nCols = length(ColNames)
  TotalCounts = matrix(0,nrow=nRows,ncol=nCols-1)

  for (CountsFile in CountsFiles) {
    writeLines(c("Adding ",CountsFile," to total counts"))

    ## Get the genotype counts for each probeset (converted or otherwise)
    if (!file.exists(CountsFile)) { stop(paste("File doesn't exist: ",CountsFile)) }
    Counts = tbl_df(read.table(CountsFile,header=TRUE,stringsAsFactors=FALSE))

    ## There should be an easier way to aggregate the counts
    if ((ncol(Counts) != nCols) || (sum(colnames(Counts) != ColNames))) {
      stop(paste("File ",CountsFile," : colnames(Counts) != ColNames"))
    }
    if ((nrow(Counts) != nRows) || (sum(Counts[[ProbeSetID]] != RowNames))) {
      stop(paste("File ",CountsFile," : rownames(Counts) != RowNames"))
    }

    Counts = select(Counts, - matches(ProbeSetID))
    Counts[is.na(Counts)] = 0
    TotalCounts = TotalCounts + Counts
  }

  TotalCounts = data.frame(RowNames,TotalCounts,stringsAsFactors=FALSE)
  colnames(TotalCounts) = ColNames
  TotalCounts = format(TotalCounts,scientific=TRUE,digits=9)
  
  Connection = file(TotalsFile,"wt")
  write.table(TotalCounts,Connection,quote=FALSE,row.names=FALSE,col.names=TRUE)
  close(Connection)
}
BatchEffectsFromCountsPerSNP <- function(ProbeSet,ps2snp,CountsFile,TotalsFile,PvalsFile,includedInTotal) {
  writeLines("BatchEffectsFromCountsPerSNP : ")
  ProbeSetID = "AffySNPID"
  
  Clusters = c("AA","AB","BB","NC")
  RowNames = ps2snp[[ProbeSetID]]
  ColNames = c(ProbeSetID,
    paste(Clusters,"all_female",sep="_"),paste(Clusters,"all_male",sep="_"),
    paste(Clusters,"ceu_female",sep="_"),paste(Clusters,"ceu_male",sep="_"),
    paste(Clusters,   "HG00097",sep="_"),paste(Clusters, "HG00264",sep="_"))
  Columns.WF = grep("ceu_female$",ColNames)
  Columns.WM = grep("ceu_male$",ColNames)
  nCols = length(ColNames)
  nRows = length(RowNames)
  Chrom = ps2snp$Chromosome
  
  if (!file.exists(TotalsFile)) { stop(paste("File doesn't exist: ",TotalsFile)) }
  if (!file.exists(CountsFile)) { stop(paste("File doesn't exist: ",CountsFile)) }

  TotalCounts = read.table(TotalsFile,header=TRUE,stringsAsFactors=FALSE)
  if ((ncol(TotalCounts) != nCols) || (sum(colnames(TotalCounts) != ColNames))) {
    stop(paste("Totals : colnames(TotalCounts) != ColNames"))
  }
  if ((nrow(TotalCounts) != nRows) || (sum(TotalCounts[[ProbeSetID]] != RowNames))) {
    stop(paste("Totals : rownames(TotalCounts) != RowNames"))
  }
  TotalCounts.WF = TotalCounts[,Columns.WF]
  TotalCounts.WM = TotalCounts[,Columns.WM]
  
  BatchCounts = read.table(CountsFile,header=TRUE,stringsAsFactors=FALSE)
  if ((ncol(BatchCounts) != nCols) || (sum(colnames(BatchCounts) != ColNames))) {
    stop(paste("Batch ",Batch," : colnames(BatchCounts) != ColNames"))
  }
  if ((nrow(BatchCounts) != nRows) || (sum(BatchCounts[[ProbeSetID]] != RowNames))) {
    stop(paste("Batch ",Batch," : rownames(BatchCounts) != RowNames"))
  }
  BatchCounts.WF = BatchCounts[,Columns.WF]
  BatchCounts.WM = BatchCounts[,Columns.WM]
  
  if (includedInTotal) {
    Statistics = apply(cbind(Chrom,TotalCounts.WF-BatchCounts.WF,TotalCounts.WM-BatchCounts.WM,
      BatchCounts.WF,BatchCounts.WM),1,CompareGenotypes)
  } else {
    Statistics = apply(cbind(Chrom,TotalCounts.WF,TotalCounts.WM,
      BatchCounts.WF,BatchCounts.WM),1,CompareGenotypes)
  }

  Statistics = data.frame(RowNames,t(Statistics))
  Statistics = format(Statistics,scientific=TRUE,digits=9)
  colnames(Statistics) = c(ProbeSetID,"pFRQ_ceu","pFRQ_ceu_female","pFRQ_ceu_male")
  write.table(Statistics,file=PvalsFile,quote=FALSE,col.names=TRUE,row.names=FALSE)
}
Males2FemalesFromCountsPerSNP <- function(ProbeSet,ps2snp,CountsFile,PvalsFile) {
  ## A comparison between males and females (genotype or allele frequencies)
  ## is very similar to the batch effects test (but with includedInTotal = FALSE)
  writeLines("Males2FemalesFromCountsPerSNP : ")
  ProbeSetID = "AffySNPID"
  
  Clusters = c("AA","AB","BB","NC")
  RowNames = ps2snp[[ProbeSetID]]
  ColNames = c(ProbeSetID,
    paste(Clusters,"all_female",sep="_"),paste(Clusters,"all_male",sep="_"),
    paste(Clusters,"ceu_female",sep="_"),paste(Clusters,"ceu_male",sep="_"),
    paste(Clusters,   "HG00097",sep="_"),paste(Clusters, "HG00264",sep="_"))
  Columns.WF = grep("ceu_female$",ColNames)
  Columns.WM = grep("ceu_male$",ColNames)
  nCols = length(ColNames)
  nRows = length(RowNames)
  Chrom = ps2snp$Chromosome
  
  if (!file.exists(CountsFile)) { stop(paste("File doesn't exist: ",CountsFile)) }

  BatchCounts = read.table(CountsFile,header=TRUE,stringsAsFactors=FALSE)
  if ((ncol(BatchCounts) != nCols) || (sum(colnames(BatchCounts) != ColNames))) {
    stop(paste("Batch ",Batch," : colnames(BatchCounts) != ColNames"))
  }
  if ((nrow(BatchCounts) != nRows) || (sum(BatchCounts[[ProbeSetID]] != RowNames))) {
    stop(paste("Batch ",Batch," : rownames(BatchCounts) != RowNames"))
  }
  BatchCounts.WF = BatchCounts[,Columns.WF]
  BatchCounts.WM = BatchCounts[,Columns.WM]
  Zeros = matrix(0,nrow=nRows,ncol=4)
  print(head(Zeros))

  ## We want to compare (BatchCounts.WF + Zeros) vs (BatchCounts.WM + Zeros)
  ## CompareGenotypes returns three p-values but only the first of those makes sense here
  Statistics = apply(cbind(Chrom,BatchCounts.WF,Zeros,Zeros,BatchCounts.WM),1,CompareGenotypes)

  Statistics = data.frame(RowNames,t(Statistics))
  ##Statistics = format(Statistics[,1:2],scientific=TRUE,digits=9)
  Statistics = format(Statistics,scientific=TRUE,digits=9)
  colnames(Statistics) = c(ProbeSetID,"pFRQ_ceu","pFRQ_ceu_female","pFRQ_ceu_male")
  write.table(Statistics,file=PvalsFile,quote=FALSE,col.names=TRUE,row.names=FALSE)
}
HardyWeinbergFromCountsPerSNP <- function(ProbeSet,ps2snp,CountsFile,PvalsFile) {
  writeLines("HardyWeinbergFromCountsPerSNP : ")
  ProbeSetID = "AffySNPID"
  
  Clusters = c("AA","AB","BB","NC")
  RowNames = ps2snp[[ProbeSetID]]
  ColNames = c(ProbeSetID,
    paste(Clusters,"all_female",sep="_"),paste(Clusters,"all_male",sep="_"),
    paste(Clusters,"ceu_female",sep="_"),paste(Clusters,"ceu_male",sep="_"),
    paste(Clusters,   "HG00097",sep="_"),paste(Clusters, "HG00264",sep="_"))
  Columns.WF = grep("ceu_female$",ColNames)
  Columns.WM = grep("ceu_male$",ColNames)
  nCols = length(ColNames)
  nRows = length(RowNames)
  Chrom = ps2snp$Chromosome
  
  if (!file.exists(CountsFile)) { stop(paste("File doesn't exist: ",CountsFile)) }

  BatchCounts = read.table(CountsFile,header=TRUE,stringsAsFactors=FALSE)
  if ((ncol(BatchCounts) != nCols) || (sum(colnames(BatchCounts) != ColNames))) {
    stop(paste("Batch ",Batch," : colnames(BatchCounts) != ColNames"))
  }
  if ((nrow(BatchCounts) != nRows) || (sum(BatchCounts[[ProbeSetID]] != RowNames))) {
    stop(paste("Batch ",Batch," : rownames(BatchCounts) != RowNames"))
  }
  BatchCounts.WF = BatchCounts[,Columns.WF]
  BatchCounts.WM = BatchCounts[,Columns.WM]
  
  Statistics = apply(cbind(Chrom,BatchCounts.WF,BatchCounts.WM),1,HardyWeinbergPvals)
  Statistics = data.frame(RowNames,t(Statistics))
  Statistics = format(Statistics,scientific=TRUE,digits=9)
  colnames(Statistics) = c(ProbeSetID,"pHWE_ceu","pHWE_ceu_female","pHWE_ceu_male")
  write.table(Statistics,file=PvalsFile,quote=FALSE,col.names=TRUE,row.names=FALSE)
}
PlateEffectsFromAxiomCallsFile <- function(ProbeSet,ps2snp,
					   BinDatadir,BatchInfo,
					   PsPerformance,PvaluesFile,
					   exclude.indivs=character(),
					   minPlateSize = 48) {
  writeLines("PlateEffectsFromAxiomCallsFile :")
  if (sub.is.ukbiobank(Batch)) {
    sub.ps2snp = ukbiobank.ps2snp(ProbeSet)
  } else if (sub.is.ukbileve(Batch)) {
    sub.ps2snp = ukbileve.ps2snp(ProbeSet)
  } else {
    stop(paste("Batch ",Batch," not found"))
  }
  
  ProbeSetID = sub.ps2snp$ProbeSetID
  Chrom = sub.ps2snp$Chromosome

  IDsToExclude = match(exclude.indivs,BatchInfo$PIID)

  CallsConnection = file(paste(BinDatadir,"/AxiomGT1.",ProbeSet,".calls.bin",sep=""),"rb")
  values = readBin(CallsConnection,numeric(),n=3)
  nProbes = as.integer(values[1])
  nSamples = as.integer(values[2])
  recSize = as.integer(values[3])
  
  writeLines(c("ProbeSet = ",ProbeSet))
  writeLines(c("nSamples = ",nSamples))
  writeLines(c("nProbesets = ",nProbes))
  fromProbe = 1
  toProbe = nProbes
  
  Gender = BatchInfo$Inferred.Gender
  if (length(Gender)!=nSamples) { stop("length(Gender)!=nSample") }
  if (!is.factor(Gender)) { Gender = as.factor(Gender) }
  
  Plate = BatchInfo$Processed.Plate.Name
  if (length(Plate)!=nSamples) { stop("length(Plate)!=nSample") }
  if (!is.factor(Plate)) { Plate = as.factor(Plate) }
  
  ## Are there plates that are less than half full?
  ## Exclude small plates because they are more likely
  ## to have different allele distribution by chance
  PlateSizes = table(Plate)
  SmallPlates = names(PlateSizes)[PlateSizes<minPlateSize]
  InSmallPlate = which(Plate %in% SmallPlates)

  ## Read the data in chunks of size = chunkSize
  probesToRead = toProbe - fromProbe + 1
  chunkSize = 10000
  firstChunk = TRUE
  
  probesToRead = toProbe - fromProbe + 1
  chunkSize = 10000
  firstChunk = TRUE
  
  U = "HG00097"
  V = "HG00264"
  uReferences = get.unique.references(BatchInfo,U)
  vReferences = get.unique.references(BatchInfo,V)
  allReferences = get.batch.references(BatchInfo)
  IncludeIndivs = 1:nSamples
  IncludeIndivs = setdiff(IncludeIndivs,allReferences)            ## always exclude the reference samples
  IncludeIndivs.AF = intersect(IncludeIndivs,which(Gender=='F'))  ## all_female
  IncludeIndivs.AM = intersect(IncludeIndivs,which(Gender=='M'))  ## all_male
  IncludeIndivs = setdiff(IncludeIndivs,IDsToExclude)             ## exclude outliers
  IncludeIndivs = setdiff(IncludeIndivs,InSmallPlate)             ## exclude "small" plates
  IncludeIndivs.WF = intersect(IncludeIndivs,which(Gender=='F'))  ## british_white_female
  IncludeIndivs.WM = intersect(IncludeIndivs,which(Gender=='M'))  ## british_white_male
  
  nFemales = length(IncludeIndivs.WF)
  nMales = length(IncludeIndivs.WM)
  writeLines(c("Number of females = ",nFemales))
  writeLines(c("Number of males = ",nMales))

  Plate = Plate[IncludeIndivs]
  Plate = as.character(Plate)
  Gender = Gender[IncludeIndivs]
  Gender = as.character(Gender)
  
  if (fromProbe>1) {
    seek(CallsConnection,(fromProbe-1)*(1*nSamples+pidSize),origin="current")
  }
  
  Pvalues = data.frame()
  writeLines("Processing...")
  while (probesToRead>0) {
    
    chunkSize = min(chunkSize,probesToRead)        
    chunkChr = Chrom[1:chunkSize] 
    Chrom = Chrom[-(1:chunkSize)]
    
    Bytes = readBin(CallsConnection,what="raw",n=chunkSize*(1*nSamples+pidSize))
    CharBytes = rep(1:pidSize,times=chunkSize) + rep(0:(chunkSize-1),each=pidSize)*(1*nSamples+pidSize)
    ProbeSets = readBin(Bytes[CharBytes],what="character",n=chunkSize)
    Calls = readBin(Bytes[-CharBytes],what="int",size=1,n=chunkSize*nSamples)
    Calls = matrix(Calls,ncol=nSamples,byrow=TRUE)
    Calls = Calls[,IncludeIndivs]
    
    chunkProbeSetID = ProbeSetID[1:chunkSize]
    ProbeSetID = ProbeSetID[-(1:chunkSize)]
    ProbeSets = sub(" *","",ProbeSets)
    
    if(sum(ProbeSets!=chunkProbeSetID)) {
      stop("sum(ProbeSets!=chunkProbeSetID)")
    }
    if (firstChunk) {
      writeLines("/* Small chunk of the calls matrix */")
      print(Calls[1:4,1:8])
    }
    
    Chunk = t(apply(cbind(chunkChr,Calls),1,calls.FisherExact,gender=Gender,category=Plate))
    Pvalues = rbind(Pvalues,data.frame(ProbeSetID=ProbeSets,Chunk,stringsAsFactors = FALSE))
    
    firstChunk = FALSE
    probesToRead = probesToRead - chunkSize
    print(toProbe - probesToRead,quote=FALSE)
  }
  
  close(CallsConnection)

  ApplyAffymetrixQCtoAffySNPs(ps2snp,PsPerformance,Pvalues,PvaluesFile)
}
