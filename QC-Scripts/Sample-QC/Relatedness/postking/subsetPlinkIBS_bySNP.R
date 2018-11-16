args = commandArgs(trailingOnly=TRUE)

library(stringr)

# args = c("/well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/Relatedness/b1__b11-b001__b095-pair_batches.filtered-duplicates-genotypes.traw","/well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/Sample-QC/Relatedness/b1__b11-b001__b095-pair_batches.filtered-duplicates-twins.txt")

genotypesFile = args[1]
dupeFile = args[2]

print("Reading in data from: ")
print(args)

outFile = gsub(".txt",".columnsInGenotypeFile.txt",dupeFile)
print(outFile)
if(outFile==dupeFile) {
    print("AAAGGH!")
    quit()
}

genotypesHeader = read.table(genotypesFile,header=FALSE,stringsAsFactors=FALSE,nrow=1)[1,]

sampleIDs = sapply(genotypesHeader[-c(1:6)],function(x){
    a = str_split(x,"_")[[1]][-c(1:3)]
    paste0(a,collapse="_")
})

dupes = read.table(dupeFile,header=TRUE,stringsAsFactors=FALSE)

dupPairs1 = paste0(dupes[,"ID1"],dupes[,"ID2"])
dupPairs2 = paste0(dupes[,"ID2"],dupes[,"ID1"])
print(head(dupPairs1))


check = sum( ! sampleIDs %in% c(dupes$ID1,dupes$ID2))
if(check>0) print("ERROR: not all your genotype samples are in the duplicates file.")
check = sum( ! c(dupes$ID1,dupes$ID2) %in% sampleIDs) 
if(check>0) print("ERROR: not all your duplicated samples are in the genotypes file.")



theseColumns = sapply( 1:nrow(dupes), function(i){
    
    id1 = dupes[i,"ID1"]    
    id2 = dupes[i,"ID2"]
    columns = c( which(sampleIDs==id1),which(sampleIDs==id2) ) + 6
#    print(cbind(c(id1,id2),c(genotypesHeader[columns])))

#    cols = myColClasses
#    cols[columns] = "numeric"
#    forSystem  = paste0("awk '{print $",columns[1],",$",columns[2],"}' ",genotypesFile," > test")
    
    #forSystem  = paste0("awk '{ if($",columns[1],"==\"NA\" || $",columns[2],"==\"NA\" ) print NR, \"NA\"; else if( $",columns[1],"!=$",columns[2]," ) print NR, \"1\" }' ",genotypesFile, " > ",dirname(genotypesFile),"/duplicateDiffs/genoDiffsDuplicates.",id1,"..",id2,".txt") # 13:57:16 to 13:56:56
    
#    print(date()) # 13:20:06 to 13:21:47
#    o = read.table(genotypesFile,colClasses = cols,header=TRUE,stringsAsFactors=FALSE)
    #print(date())
    #system(forSystem,wait=TRUE) # 13:27:01 to 13:27:23 = 22 secs
#   # o2 = read.table('test',colClasses="integer",header=TRUE,stringsAsFactors=FALSE)
    #print(date())
    return(columns)
})

toPrint = cbind(dupes[,c("ID1","ID2")],t(theseColumns))

write.table(toPrint,file=outFile,quote=FALSE,col.names=FALSE,row.names=FALSE)

print(paste0("List of columns written to: ",outFile))
