

source("/well/ukbiobank/qcoutput.V2_QCed.sample-QC-testing/QC-Scripts/R/scripts/ukbbcolors.R")


sub.is.ukbileve <- function(Batches) {
  x = rep(FALSE,length(Batches))
  x[grep("UKBiLEVEAX_b",Batches)] = TRUE
  return(x)
}
sub.is.ukbiobank <- function(Batches) {
  x = rep(FALSE,length(Batches))
  x[grep("Batch_b",Batches)] = TRUE
  return(x)
}
sub.plate.info <- function(Batch,Info) {
  ## It used to be that the sample tables from Affymetrix contained only one plate name
  ## (so I did not know it corresponds to the submitted plate name)
  if (sub.is.ukbileve(Batch)) {
    Plate = Info$Submitted.Plate.Name
  } else if (sub.is.ukbiobank(Batch)) {
    Plate = Info$Processed.Plate.Name
  }
  return(Plate)
}
sub.well.info <- function(Batch,Info) {
  if (sub.is.ukbileve(Batch)) {
    Well = Info$Submitted.Well
  } else if (sub.is.ukbiobank(Batch)) {
    Well = Info$Processed.Well
  }
  return(Well)
}

## OLD - from interim release
#sub.assign.sample.piid <- function(Batch,Info) {
#  Sample = Info$Sample.Name
#  Plate = sub.plate.info(Batch,Info)
#  Well = sub.well.info(Batch,Info)
#  if (Batch=="UKBiLEVEAX_b11") {
#    PIID = paste(Plate,Sample,Well,sep="__")
#  } else {
#    PIID = paste(Plate,Sample,sep="__")
#  }
#  return(PIID)
#}

sub.assign.sample.piid <- function(Batch,Info) {
    PIID = Info$Best.Array
    return(PIID)
}

    
sub.gender.factor <- function(Gender) {
  Gender = as.character(Gender)
  is.male = which((Gender=="Male")|(Gender=="male"))
  is.female = which((Gender=="Female")|(Gender=="female"))
  n = length(Gender)
  Gender = rep("0",n)
  Gender[is.male] = "M"
  Gender[is.female] = "F"
  Gender = as.factor(Gender)
  return(Gender)
}
get.batch.no <- function(Batch) {
  m = regexec("_(b[0-9]+)",Batch)
  t = regmatches(Batch,m)[[1]]
  x = t[2]
  return(x)
}
get.batch.info <- function(Batch,CsvFile) {

  Info = read.csv(CsvFile,stringsAsFactors=FALSE)

  ## Exclude samples that have failed Affymetrix sample QC
  Info = Info[Info$Sample.Status=="Pass",]
  ## Order alphabetically by Best Array filename
  ## This is the order of calls/confidences/intensities in the binary files
  Info = Info[order(Info$Best.Array),]
  
  Pops = ethnicity2pop(Info$Ethnic.background)
  Chars = ethnicity2char[Pops]
  Colors = ethnicity2col[Pops]
  Info$PIID = sub.assign.sample.piid(Batch,Info)
  Info$Batch = Batch
  Info$Pops = Pops
  Info$Chars = Chars
  Info$Colors = Colors
  Info$Submitted.Gender = sub.gender.factor(Info$Submitted.Gender)
  Info$Inferred.Gender = sub.gender.factor(Info$Inferred.Gender)

  return(Info)
}
#get.unique.references <- function(Info,x) {
#  Sample = Info$PIID
#  References = grep(x,Sample)
#  return(References)
#}

# Clare's update
get.unique.references <- function(Info,x) {
  Sample = Info$Sample.Name
  References = grep(x,Sample)
  return(References)
}


get.batch.references <- function(Info) {
  References = which(Info$Sample.Source!="Customer")
  return(References)
}

### clare's additions (from mybatchinfo.R)
ukbileve.batches <- function( ) {
  Batches = paste("UKBiLEVEAX_b",1:11,sep="")
  return(Batches)
}
ukbiobank.batches <- function( ) {
  ## It seems that I can't get my perl script to work on the summary file for b072
  ## There is an error about an invalid archive.
  #Batches = paste("Batch_b",sprintf("%03d",c(1:71,73:80)),sep="")
   Batches = paste("Batch_b",sprintf("%03d",c(1:94)),sep="")
  return(Batches)
}
ukbiobank.batches.release1 <- function( ) {
  Batches = paste("Batch_b",sprintf("%03d",1:22),sep="")
  return(Batches)
}
all.batches <- function() {
    Batches = c( ukbileve.batches(),ukbiobank.batches() )
    return(Batches)
}
