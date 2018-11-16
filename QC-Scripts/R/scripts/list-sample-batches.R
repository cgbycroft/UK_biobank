#######
# script to create a file that lists batch for each sample (effectively just cat .fam files together)
#######

library(stringr)

args = commandArgs(trailingOnly=T)

h = args[-c(which(args%in%c("-out")),1+which(args%in%c("-out")))]
for(helperScript in h){
    source(helperScript)
}

outfile = args[which(args=="-out")+1]

b=1
for( Batch in all.batches() ){

    print(Batch)
    famFile = paste0(baseDataDir,"/",Batch,"/sorted_bed/",Batch,"-sexchrom.fam")
    fam = read.table(famFile,stringsAsFactors=FALSE,colClasses=c("character",rep("NULL",5)) )
    fam$batch = Batch
    
    if(b==1) write.table(fam,file=outfile,quote=FALSE,row.names=FALSE,col.names=FALSE) else write.table(fam,file=outfile,quote=FALSE,row.names=FALSE,col.names=FALSE,append=TRUE)

    b=2
}

print(paste0("DONE. See ",outfile))
