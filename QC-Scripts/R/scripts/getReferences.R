#######
# script to create file with list of reference samples for sample QC, with BestArray; as well as a file that lists batch for each sample (just cat .fam files together)
#######

library(stringr)

args = commandArgs(trailingOnly=T)

h = args[-c(which(args%in%c("-out")),1+which(args%in%c("-out")))]
for(helperScript in h){
    source(helperScript)
}

outfile = args[which(args=="-out")+1]

hell = unlist(list.all.references())  # this function is written in auxFunctions.R

print(length(hell))

write.table(hell,file=outfile,quote=FALSE,row.names=FALSE,col.names=FALSE)
