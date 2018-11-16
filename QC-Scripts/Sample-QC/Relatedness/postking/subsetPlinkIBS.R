args = commandArgs(trailingOnly=TRUE)


ibdFile = args[1]
dupeFile = args[2]

print("Reading in data from: ")
print(args)

ibd = read.table(ibdFile,header=TRUE,stringsAsFactors=FALSE)
dupes = read.table(dupeFile,header=TRUE,stringsAsFactors=FALSE)

IbdPairs = paste0(ibd[,2],ibd[,4])
print(head(IbdPairs))

dupPairs1 = paste0(dupes[,"ID1"],dupes[,"ID2"])
dupPairs2 = paste0(dupes[,"ID2"],dupes[,"ID1"])
print(head(dupPairs1))

newIbd = ibd[IbdPairs%in%c(dupPairs1,dupPairs2),]


if( dim(newIbd)[1]!=dim(dupes)[1] ) print("number of pairs in new ibd file doesn't match number of pairs in the dupes file... Check output.")

write.table(newIbd,file=paste0(ibdFile,".dupes"),quote=FALSE,col.names=TRUE,row.names=FALSE)

print(paste0("Subset file written to: ",ibdFile,".dupes"))
