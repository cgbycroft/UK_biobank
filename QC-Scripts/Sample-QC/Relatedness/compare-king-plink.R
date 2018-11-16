## Script to compare results from plink and king (i.e the overlap)


args = commandArgs(trailingOnly=T)
#args = c("/well/ukbiobank/expt/V2_QCed.SNP-QC/src/V2_QCed.snpqc-tests.R","/well/ukbiobank/expt/V2_QCed.SNP-QC/src/V2_QCed.bin2clusterplots.R","/well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/R/scripts/readPSperformance.R","/well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/R/scripts/auxFunctions.R","-ibd" ,"/well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/Relatedness/b1__b11-b001__b095-autosome-sampleqc-filtered.0.003-White_British.genome.combined","-rel","/well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/Relatedness/b1__b11-b001__b095-pair_batches.filtered-samples-pruned-200-White_British.kin0","-ibs0","0.0012")

print(args)
h = args[-c(which(args%in%c("-rel","-ibd","-ibs0","-hm")),1+which(args%in%c("-rel","-ibd","-ibs0","-hm")))]
for(helperScript in h){
    source(helperScript)
}

PlinkOutFile = args[which(args=="-ibd")+1]
KingOutFile = args[which(args=="-rel")+1]
OutputFile = "compare-king-plink"
ibs0Threshold = as.numeric(args[which(args=="-ibs0")+1])

#het missing information
load(paste0(baseSampleQCDir,"/data/Combined/b1__b11-b001__b095-autosome-sampleqc-hetcorrected-imiss.RData"),verbose=TRUE)

# read kinship files
kinPlink = read.table(PlinkOutFile,header=TRUE,stringsAsFactors=FALSE)  # this is a huge file...
kinKing = read.table(KingOutFile,header=TRUE,stringsAsFactors=FALSE)

# get classes
classesKing = get.kin.classes(kinKing,ibs0Threshold=ibs0Threshold)
kinKing$class = classesKing
classesPlink = get.kin.classes(kinPlink,ibs0Threshold=ibs0Threshold)
kinPlink$class = classesPlink

print("Finding overlapping pairs...")
                                        # find overlapping pairs
hm = left_join(kinKing,kinPlink,by=c("ID1"="ID1","ID2"="ID2"))
hm2 = left_join(kinKing,kinPlink,by=c("ID1"="ID2","ID2"="ID1"))

kinOverlap = rbind(hm[!is.na(hm$FID1.y),],hm2[!is.na(hm2$FID1.y),])
colnames(kinOverlap)[1:9] = colnames(kinKing)[1:9]

kinOverlapPlink = kinOverlap[,c(1:4,12,13,7,15)] # use same IBS0 as king
kinOverlapPlink$Kinship = kinOverlapPlink$Kinship.y
classesPlink2 = get.kin.classes(kinOverlapPlink,ibs0Threshold=ibs0Threshold) # where use the same IBS0 as king
kinOverlap$class.y2 = classesPlink2  # ====> no difference with class.y

print( paste0(dim(kinOverlap)[1]," pairs in overlap"))
print( paste0( "= ",100*dim(kinOverlap)[1]/dim(kinKing)[1],"% of king results."))


                                        # breakdown by kinship class
print("King all")
print( table(kinKing$class) )# king all
print("Plink all")
print( table(kinPlink$class) )# plink all

print("King overlap")
table(kinOverlap$class) # king overlap
print("Plink overlap")
table(kinOverlap$class.y) # plink overlap

# mismatches
mismatches = kinOverlap[kinOverlap$class!=kinOverlap$class.y,]

kinKingNotInOverlap = kinKing[!paste0(kinKing$ID1,kinKing$ID2)%in%paste0(kinOverlap$ID1,kinOverlap$ID2),]

print("king not in overlap")
print( table(kinKingNotInOverlap$class) ) 

outFile=paste0(gsub(".combined","",basename(PlinkOutFile)),"-WITH-",gsub(".kin0","",basename(KingOutFile)))
write.table( kinOverlap,file= paste0(outFile,".overlap.txt"))
write.table( kinKingNotInOverlap,file= paste0(outFile,".kingNotPlink.txt") )

quit()


######## PLOT

colors = colDef[kinOverlap$class]
ord = order.by.number.occurrences(colors)



png(paste0("plots/",OutputFile,"-hists.png"),width=1000,height=1000,res=150)
par(mfrow=c(2,1))
hist(kinOverlap$Kinship,breaks=1000,xlab="Kinship - KING",xlim=c(0,0.5))
abline(v=c(1/2,1/4,1/8,1/16),col="red")
hist(kinOverlap$Kinship.y,breaks=1000,xlab="Kinship - Plink",xlim=c(0,0.5))
abline(v=c(1/2,1/4,1/8,1/16),col="red")
dev.off()

                                        # scatter plot for kinship, with colours based on king
x = kinOverlap$Kinship
y = kinOverlap$Kinship.y

png(paste0("plots/",OutputFile,"-Kinship.png"),width=1000,height=1000,res=150)
plot(x[ord],y[ord],xlab="King",ylab="Plink",col=add.alpha(colors[ord],0.5),pch=16,xlim=range(c(x,y)),ylim=range(c(x,y)))
abline(0,1,col="black")
legend("bottomright",legend=names(colDef)[c(2,5,1,3,4)],col=colDef[c(2,5,1,3,4)],pch=16,bty="n")
dev.off()



                                        # relationship with missing rates??
y = (kinOverlap$Kinship.y - kinOverlap$Kinship)/kinOverlap$Kinship
x = kinOverlap$N_SNP/max(kinOverlap$N_SNP) 
    
png(paste0("plots/",OutputFile,"-KinshipDifferenceByN_SNPs.png"),width=1000,height=1000,res=150)
plot(x[ord],y[ord],xlab="fraction missing SNPs",ylab="Kinship estimate difference (Plink - KING)",col=add.alpha(colors[ord],0.5),pch=16)
legend("bottomright",legend=names(colDef)[c(2,5,1,3,4)],col=colDef[c(2,5,1,3,4)],pch=16,bty="n")
dev.off()




######### Investigate weird cluster
                                        # isolate cluster of 3rd-degree relatives. What's going on with them??
weird = ((kinOverlap$Kinship.y - kinOverlap$Kinship)/kinOverlap$Kinship) > 0.5
WeirdKin = kinOverlap[weird,]

                                        #add in het information
index1 = match(kinOverlap$ID1,Table$IID)
index2 = match(kinOverlap$ID2,Table$IID)

# heterozygosity difference? <==== Nope!
x = abs(Table$het.corrected[index1] - Table$het.corrected[index2])
y = (kinOverlap$Kinship.y - kinOverlap$Kinship)/kinOverlap$Kinship

png(paste0("plots/",OutputFile,"-KinshipDifferenceByHetDifference.png"),width=1000,height=1000,res=150)
plot(x[ord],y[ord],xlab="Heterozygosity difference (all SNPs)",ylab="Kinship estimate difference (Plink - KING)",col=add.alpha(colors[ord],0.5),pch=16)
#legend("bottomright",legend=names(colDef)[c(2,5,1,3,4)],col=colDef[c(2,5,1,3,4)],pch=16,bty="n")
dev.off()
