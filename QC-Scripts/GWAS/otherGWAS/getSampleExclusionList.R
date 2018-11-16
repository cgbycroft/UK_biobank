args <- commandArgs(trailing=TRUE)

# Read in .fam files, or otherwise. As long as the first two columns are Sample IDs.
A = read.table(args[1],header=FALSE,stringsAsFactors=FALSE)
B = read.table(args[2],header=FALSE,stringsAsFactors=FALSE)

toExclude = A[!A[,1]%in%B[,1],1:2] # Who is in A but not in B?

print(paste0(nrow(toExclude)," samples in A but not in B."))

write.table(toExclude,file=args[3],quote=FALSE,col.names=FALSE,row.names=FALSE)
