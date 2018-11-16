 # This script takes the output from ./postking/duplicate-concordance.sh and plots summaries of the concordance rates

args = commandArgs(TRUE)

#args = c("/well/ukbiobank/expt/V2_QCed.SNP-QC/src/V2_QCed.snpqc-tests.R","/well/ukbiobank/expt/V2_QCed.SNP-QC/src/V2_QCed.bin2clusterplots.R","/well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/R/scripts/readPSperformance.R","/well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/R/scripts/auxFunctions.R","-in","b1__b11-b001__b095-pair_batches.filtered-samples-pruned-200")

print(args)
h = args[-c(which(args%in%c("-in")),1+which(args%in%c("-in")))]
for(helperScript in h){
    source(helperScript)
}

setwd(paste0(baseSampleQCDir,"/QC-Scripts/Sample-QC/Relatedness"))

inPrefix = args[which(args=="-in")+1]

# raw kinship data
#kinRaw = read.table('../../../data/Relatedness/b1__b11-b001__b095-pair_batches.filtered.kin0',header=TRUE,stringsAsFactors=FALSE)

# Raw sample tables (before Oxford QC)
otherInfo = read.multiple.batch.info(batchInfoFields)


# All the duplicates and twins
twinDupeFile = paste0(inPrefix,"-duplicates-twins.txt")
twinDupe = read.table(twinDupeFile,header=TRUE,stringsAsFactors=FALSE)

# Annotated duplicates list. This was based on the 894 pairs detected in the pre-pruned kinship table (as in kinRaw).
dupeAnnotFile = "160513-WCSG-duplicated-samples.txt"
dupeAnnot = read.delim(dupeAnnotFile,header=TRUE,stringsAsFactors=FALSE,sep="\t")


# Master blind-spiked duplicates list from Sam Welsh (nee Murphy)
dupeBlindFile = "170221-BSD-pairs-DISTRIBUTED-v2.csv"
dupeBlind = read.delim(dupeBlindFile,header=TRUE,stringsAsFactors=FALSE,sep=",")
dupeBlind$Best.Array.caps = toupper(dupeBlind$Best.Array)
dupeBlind$Sibling.Best.Array.caps = toupper(dupeBlind$Sibling.Best.Array)

sib1 = dupeBlind[,!grepl("Sibling",colnames(dupeBlind))][,-c(5,8)]
sib2 = dupeBlind[,grepl("Sibling",colnames(dupeBlind))]


# Exclude any pairs with a sample that isn't in our sample tables.
dupeBlind2 = dupeBlind[(sib1[,7]%in%otherInfo$PIID)&(sib2[,7]%in%otherInfo$PIID),]
print("Pairs of intended duplicates where both best.array ids are actually in our raw sample tables:")
print(dim(dupeBlind2)[1])

# the reverse (what's left?)
dupeBlind3 =  dupeBlind[!((sib1[,7]%in%otherInfo$PIID)&(sib2[,7]%in%otherInfo$PIID)),]

j = dupeAnnot[intended,]
inAnnot = dupeBlind3[(dupeBlind3$Sibling.Best.Array.caps%in%c(dupeAnnot$Best_Array.x,dupeAnnot$Best_Array.y))|(dupeBlind3$Best.Array.caps%in%c(dupeAnnot$Best_Array.x,dupeAnnot$Best_Array.y)),]

inAnnot$myDuplicate = twinDupe$ID2[match(inAnnot$Sibling.Best.Array.caps,twinDupe$ID1)]
inAnnot$myDuplicate[is.na(inAnnot$myDuplicate)] = twinDupe$ID1[match(inAnnot$Sibling.Best.Array.caps[is.na(inAnnot$myDuplicate)],twinDupe$ID2)]
inAnnot$myDuplicate[is.na(inAnnot$myDuplicate)] = twinDupe$ID2[match(inAnnot$Best.Array.caps[is.na(inAnnot$myDuplicate)],twinDupe$ID1)]
inAnnot$myDuplicate[is.na(inAnnot$myDuplicate)] = twinDupe$ID1[match(inAnnot$Best.Array.caps[is.na(inAnnot$myDuplicate)],twinDupe$ID2)]

# are there duplicates intended, but 


# Did we find all of them???
dupeBlind2Pairs = c(paste0(dupeBlind2$Best.Array.caps,dupeBlind2$Sibling.Best.Array.caps),paste0(dupeBlind2$Sibling.Best.Array.caps,dupeBlind2$Best.Array.caps))
twinDupeP = c(paste0(twinDupe$ID2,twinDupe$ID1),paste0(twinDupe$ID1,twinDupe$ID2))

notFound = !paste0(dupeBlind2$Best.Array.caps,dupeBlind2$Sibling.Best.Array.caps)%in%twinDupeP
missed = dupeBlind2[notFound,]
missedUnique = unique(c(missed$Sibling.Best.Array.caps,missed$Best.Array.caps))
sum(!missedUnique%in%otherInfo$PIID)

########
# Check set of pairs by Sample ID
########

sampleNameTab  = table(otherInfo$CollectionID)
dupSampleNames = names(sampleNameTab)[sampleNameTab>1]
dupSampleInfo = otherInfo[otherInfo$CollectionID%in%dupSampleNames,]
print(length(dupSampleNames))  # 585 unique names
print(dim(dupSampleInfo)[1])  # 1179 pairs

pairsFound = sapply(unique(dupSampleInfo$CollectionID),function(s){
    inds = dupSampleInfo$Best.Array[dupSampleInfo$CollectionID==s]
    ps = sum((twinDupe$ID1%in%inds)&(twinDupe$ID2%in%inds))
    any = sum((twinDupe$ID1%in%inds)|(twinDupe$ID2%in%inds))
    c(ps,any)
})

which(pairsFound[2,]==0)
# 2 duplicated sample names with no individuals in our duplicates list
# 5 duplicated sample names with only ONE individual in our duplicates list

# UKBL__2279896
# UKBL__4182533



# Have we found extras?
dupeAnnotIntended = dupeAnnot[intended,]
dupeAnnotIntended$Best.Array.x%in%dupeBlind2$Best.Array.

# Colin's duplicates exclusion list (these are excluded from the release output)
dupeExclFile = "/well/ukbiobank/expt/V2_QCed.identical_samples/data/V2_QCed.duplicates_exclude.txt"
dupeExl = read.delim(dupeExclFile,header=FALSE,stringsAsFactors=FALSE,sep="\t")[,1]

#dupeAnnot[(dupeAnnot$Best_Array.x%in%dupeExl)&(dupeAnnot$Best_Array.y%in%dupeExl),"Status..SM."]

unintended = dupeAnnot$Status..SM.=="Unintended"
intended = dupeAnnot$Status..SM.=="Intended"
unintendedPairs = c(paste0(dupeAnnot$Best_Array.x[unintended],dupeAnnot$Best_Array.y[unintended]),paste0(dupeAnnot$Best_Array.y[unintended],dupeAnnot$Best_Array.x[unintended]))
intendedPairs = c( paste0(dupeAnnot$Best_Array.x[intended],dupeAnnot$Best_Array.y[intended]),paste0(dupeAnnot$Best_Array.y[intended],dupeAnnot$Best_Array.x[intended]))
twn = dupeAnnot$Status..SM.=="Twins"
confirmedTwinPairs = c( paste0(dupeAnnot$Best_Array.x[twn],dupeAnnot$Best_Array.y[twn]),paste0(dupeAnnot$Best_Array.y[twn],dupeAnnot$Best_Array.x[twn]))


#missingDupes = dupeAnnot[(!dupeAnnot$Best_Array.y%in%twinDupe$ID1)&(!dupeAnnot$Best_Array.y%in%twinDupe$ID2),]
#nKin = table(c(kinRaw$ID1,kinRaw$ID1))
#nKin[missingDupes$Best_Array.x]
#nKin[missingDupes$Best_Array.y]

# NOTE: There are four pairs of duplicates which end up being excluded after pruning the kinship table.
# The following two are in the released data, but filtered out of kinship:
# A550484-4254432-072416-927_B09  A550484-4239209-122815-940_D04   <= twins. But *D04 was in the list of hetmiss outliers. Both individuals are in the release output data.
# A550465-4276624-031217-258_G03  A550465-4196233-091014-716_H06 <== Not twins, but A550465-4196233-091014-716_H06 has > 200 other 'relatives'.

# The other two pairs involve one individual that's excluded from the release output data. They are unintended duplicates, and have one in common: A550465-4195511-082814-525_E05. This is in Colin's duplicates exclusion list. Its 'partner' duplicates are themselves duplicates.
# Furthermore, A550465-4195511-082814-525_E05 has > 200 other 'relatives'.


# Who are genuine twins? This has to be generated by find-families.R Might not exist...
#realTwinsFile = paste0(inPrefix,"-nodupes-duplicates-twins.txt")
#realTwins = read.table(realTwinsFile,header=TRUE,stringsAsFactors=FALSE)
realTwins  = twinDupe[(!twinDupe$ID1 %in% dupeExl)&(!twinDupe$ID2 %in% dupeExl),] # exclude any pair with one individual that's in the exclusion file.
# Find the duplicate pairs
twins = paste0(twinDupe$ID1,twinDupe$ID2)%in%paste0(realTwins$ID1,realTwins$ID2) 
dupes = !twins

dim(realTwins)[1]==sum(twins)


# Read in the results from plink --genome output

ibdFile = paste0(baseSampleQCDir,"/data/Relatedness/",inPrefix,"-duplicates-genetic-distances-release.genome.dupes")
ibd = read.table(ibdFile,header=TRUE,stringsAsFactors=FALSE)

# count snps that are exactly the same genotype (IBS2)
IBS2 = sapply(1:nrow(twinDupe),function(i){
    s1 = twinDupe[i,"ID1"]
    s2 = twinDupe[i,"ID2"]
    S = which(( (ibd[,"IID1"]==s1)&(ibd[,"IID2"]==s2) ) | ( (ibd[,"IID2"]==s1)&(ibd[,"IID1"]==s2) ))
    d = ibd[S,"IBS2"]
    #print(d)
    return(d)
})

# count snps that are non-missing in both
NSNPs = sapply(1:nrow(twinDupe),function(i){
    s1 = twinDupe[i,"ID1"]
    s2 = twinDupe[i,"ID2"]
    S = which(( (ibd[,"IID1"]==s1)&(ibd[,"IID2"]==s2) ) | ( (ibd[,"IID2"]==s1)&(ibd[,"IID1"]==s2) ))
    d = sum(ibd[S,c("IBS2","IBS1","IBS0")])
    #print(d)
    return(d)
})


fracs = IBS2/NSNPs

twinDupe$IBS2 = IBS2
twinDupe$NMsnps = NSNPs
twinDupe$fracConcordance = fracs
twinDupe$fracDiscordance = 1-fracs

# Subset to just the duplicates
dupePairs = twinDupe[dupes,]
twinPairs = twinDupe[!dupes,]

# write out the results
write.table(twinDupe,file= paste0(inPrefix,"-duplicates-twins-with-discordance.txt"),quote=FALSE,col.names=TRUE,row.names=FALSE)

write.table(dupePairs,file= paste0(inPrefix,"-duplicates-with-discordance.txt"),quote=FALSE,col.names=TRUE,row.names=FALSE)

#########
save(twinDupe,dupePairs,realTwins,dupeAnnot,file=paste0(inPrefix,"-duplicates-twins-with-discordance.RData"))
#########

                                        # plot the results
print(paste0(dim(dupePairs)[1]," duplicated pairs in kinship table."))
print(paste0(sum(twins)," genuine twin pairs in kinship table."))

print(paste0(length(unique(c(dupePairs$ID1,dupePairs$ID2)))," unique samples among duplicate pairs."))
print(paste0("Mean fraction discordance: ",mean(dupePairs$fracDiscordance,na.rm=TRUE)))
print(paste0("SD fraction discordance: ",sd(dupePairs$fracDiscordance,na.rm=TRUE)))
print(paste0("Range fraction discordance: ",paste(range(dupePairs$fracDiscordance,na.rm=TRUE),collapse=" to ")))



# Duplicates (I.e not in twins after applying exclusion list)
png(paste0("plots/",inPrefix,"-duplicates-concordance.png"),width=1000,height=1000,res=150)
hist(100*dupePairs$fracConcordance,breaks=30,xlab="% of non-missing genotypes identical",col="darkgray",main=paste0("Concordance rates for ",dim(dupePairs)[1]," pairs of duplicated samples"),
cex=2,xlim=100*c(min(twinDupe$fracConcordance),1))
dev.off()

png(paste0("plots/",inPrefix,"-duplicates-discordance.png"),width=1000,height=1000,res=150)
hist(100*dupePairs$fracDiscordance,breaks=30,xlab="% of non-missing genotypes discordant",col="darkgray",main=paste0("Discordance rates for ",dim(dupePairs)[1]," pairs of duplicated samples"),
cex=2,xlim=100*c(0,max(twinDupe$fracDiscordance)))
dev.off()


# Unintended uplicates based on the annotated duplicates file
these=paste0(twinDupe$ID1,twinDupe$ID2)%in%unintendedPairs
png(paste0("plots/",inPrefix,"-duplicates-unintended-concordance.png"),width=1000,height=1000,res=150)
hist(100*twinDupe$fracConcordance[these],breaks=30,xlab="% of non-missing genotypes identical",col="darkgray",main=paste0("Concordance rates for ",dim(twinDupe[these,])[1]," pairs of unintended duplicated samples"),
cex=2,xlim=100*c(min(twinDupe$fracConcordance),1))
dev.off()

png(paste0("plots/",inPrefix,"-duplicates-unintended-discordance.png"),width=1000,height=1000,res=150)
hist(100*twinDupe$fracDiscordance[these],breaks=30,xlab="% of non-missing genotypes discordant",col="darkgray",main=paste0("Discordance rates for ",dim(twinDupe[these,])[1]," pairs of unintended duplicated samples"),
cex=2,xlim=100*c(0,max(twinDupe$fracDiscordance)))
dev.off()



# Intended uplicates based on the annotated duplicates file
these=paste0(twinDupe$ID1,twinDupe$ID2)%in%intendedPairs
png(paste0("plots/",inPrefix,"-duplicates-intended-concordance.png"),width=1000,height=1000,res=150)
hist(100*twinDupe$fracConcordance[these],breaks=30,xlab="% of non-missing genotypes identical",col="darkgray",main=paste0("Concordance rates for ",dim(twinDupe[these,])[1]," pairs of intended duplicated samples"),
cex=2,xlim=100*c(min(twinDupe$fracConcordance),1))
dev.off()

png(paste0("plots/",inPrefix,"-duplicates-intended-discordance.png"),width=1000,height=1000,res=150)
hist(100*twinDupe$fracDiscordance[these],breaks=30,xlab="% of non-missing genotypes discordant",col="darkgray",main=paste0("Discordance rates for ",dim(twinDupe[these,])[1]," pairs of intended duplicated samples"),
cex=2,xlim=100*c(0,max(twinDupe$fracDiscordance)))
dev.off()


# Twins
png(paste0("plots/",inPrefix,"-twins-concordance.png"),width=1000,height=1000,res=150)
hist(100*twinPairs$fracConcordance,breaks=30,xlab="% of non-missing genotypes identical",col="darkgray",main=paste0("Concordance rates for ",dim(twinPairs)[1]," pairs of twins"),
cex=2,xlim=100*c(min(twinDupe$fracConcordance),1))
dev.off()

png(paste0("plots/",inPrefix,"-twins-discordance.png"),width=1000,height=1000,res=150)
hist(100*twinPairs$fracDiscordance,breaks=30,xlab="% of non-missing genotypes discordant",col="darkgray",main=paste0("Discordance rates for ",dim(twinPairs)[1]," pairs of twins"),
cex=2,xlim=100*c(0,max(twinDupe$fracDiscordance)))
dev.off()


