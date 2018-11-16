#######
# script to create functions that get useful auxiliary information regarding ps performance

# PsPerformace is a string giving the ps.performance file to read

read.ps.performance <- function(PsPerformance,columns) {

  columns = c("probeset_id","affy_snp_id",columns)
  header = read.table(gzfile(PsPerformance),header=FALSE,stringsAsFactors=FALSE,nrow=1)[1,]
  colClasses = rep("NULL",length(header))
  colClasses[header%in%columns] = NA
  
  filtered = read.table(gzfile(PsPerformance),header=TRUE,stringsAsFactors=FALSE,colClasses=colClasses)
  filtered = tbl_df(filtered)

  return(filtered)
}
filter.ps.performance <- function(PsPerformance,column,value) {
    # this function returns a data.frame of SNPs (affy_snp_id) and <column>, where all snps should be <value> in the specified column
    header = read.table(gzfile(PsPerformance),header=FALSE,stringsAsFactors=FALSE,nrow=1)[1,]
    colNumber = which(header==column)
    forSystem = paste0("zcat ",PsPerformance," | awk '$",colNumber," == ",value," {print $1,$2,$",colNumber,"}' > tmp")
    system(forSystem)
    filtered = read.table("tmp",header=FALSE,stringsAsFactors=FALSE)
    system("rm tmp")
    return(filtered)
}
