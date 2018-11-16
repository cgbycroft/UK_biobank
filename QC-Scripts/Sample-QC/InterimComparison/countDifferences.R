args = commandArgs(trailingOnly=TRUE)

#args = c("/well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/Combined/interimReleaseSamples/V2_QCed.export.UKBiLEVEAX_b1-b11.Batch_b001_b095.chr1.interim_samples_UKBiLEVEAX_b1.traw","UKBiLEVEAX_b1.1.new.best.array.txt","/well/ukbiobank/expt/check_data_release_1_ctsu/data/calls/UKBL_UKBiLEVEAX_b1.chr1.calls.ctsu_basic.txt","UKBiLEVEAX_b1.1.interim.best.array.txt")

#args = c("/well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/Combined/interimReleaseSamples/V2_QCed.export.UKBiLEVEAX_b1-b11.Batch_b001_b095.chr22.interim_samples_UKBiLEVEAX_b1.traw","/well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/Combined/interimReleaseSamples/Comparison/UKBiLEVEAX_b1.22.new.best.array.txt","/well/ukbiobank/expt/check_data_release_1_ctsu/data/calls/UKBL_UKBiLEVEAX_b1.chr22.calls.ctsu_basic.txt","/well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/Combined/interimReleaseSamples/Comparison/UKBiLEVEAX_b1.22.interim.best.array.txt")

print(args)

h = c("/well/ukbiobank/expt/V2_QCed.SNP-QC/src/V2_QCed.snpqc-tests.R","/well/ukbiobank/expt/V2_QCed.SNP-QC/src/V2_QCed.bin2clusterplots.R","/well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/R/scripts/readPSperformance.R","/well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/R/scripts/auxFunctions.R")
for(s in h) source(s)

# snp annotations in new release
type="autosome"
if(grepl("chrX|chrY|chrMT",args[1]) ) type="sexchrom"
print(type)
BB.ps2snp = ukbiobank.ps2snp(type)
BL.ps2snp = ukbileve.ps2snp(type)
ps2snpAll=rbind(BB.ps2snp,BL.ps2snp)
#ps2snpAll=rbind(BB.ps2snpSex,BL.ps2snpSex)
ps2snpAll=unique(ps2snpAll)


# frequency calculations (from compute-frequencies.sh)
snpFrequencyData = "/well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/Combined/V2_QCed.export.UKBiLEVEAX_b1-b11.Batch_b001_b095.chr"

outFile = paste0(args[5],".RData")

# Align the samples
newFamRaw = read.table(args[2],header=FALSE,stringsAsFactors=FALSE)[,1] # new data sample order
underscores = 1+(str_count(newFamRaw,"_")-1)/2
sum(underscores==3)
newFam = sapply(1:length(newFamRaw),function(x) paste(str_split(newFamRaw[x],"_")[[1]][1:underscores[x]],collapse="_") )

interimFam = read.table(args[4],header=FALSE,stringsAsFactors=FALSE)[,1] # interim data sample order

sampleMatch = match(interimFam,newFam) # order of the interim samples in coordinates of new samples
sampleOrder = sampleMatch[!is.na(sampleMatch)]

intersectSamples =  intersect(interimFam,newFam)

print( paste0( length(intersectSamples)," samples in intersection." ) )
print( paste0( sum(!interimFam%in%newFam)," samples not in new data." ) )
print( paste0( sum(!newFam%in%interimFam)," samples not in old data." ) )



# Align the SNPs
newSnps = read.table(gsub("best\\.array","snp\\.order",args[2]),header=FALSE,stringsAsFactors=FALSE)[,1]
interimSnps = read.table(gsub("best\\.array","snp\\.order",args[4]),header=FALSE,stringsAsFactors=FALSE)[,1]


# convert new snps to Affy IDs
newSnps2 = gsub("_.*","",newSnps)
if( sum( !newSnps2%in%c(ps2snpAll$AffySNPID,ps2snpAll$dbSNPRSID) )!=0 ) print("AGH, can't find all the snps in the new data in the ps2snp files!")
newSnps2 = ps2snpAll$AffySNPID[match(newSnps2,ps2snpAll$dbSNPRSID)]
newSnps2[is.na(newSnps2)] = newSnps[is.na(newSnps2)]

intersectSnps = intersect(interimSnps,newSnps2)

snpMatch = match(interimSnps[interimSnps%in%intersectSnps],newSnps2) # order of the interim snps in coordinates of new snps
snpOrder = snpMatch[!is.na(snpMatch)]

snpOrderOld = interimSnps%in%intersectSnps

print( paste0( length(intersectSnps)," snps in intersection." ) )
print( paste0( sum(!interimSnps%in%newSnps2)," snps not in new data." ) )
print( paste0( sum(!newSnps2%in%interimSnps)," snps not in old data." ) )


###### Read in all of the new release data since it is in long form

print("Reading in and processing new 500K data for the intersection samples...")
nSamples = length(newFam)

colclasses = c("numeric","character","numeric","numeric","character","character",rep("numeric",nSamples))

# takes a while to read in the data for a whole batch, but saves time later
newData = read.table(args[1],header=TRUE,stringsAsFactors=FALSE,colClasses=colclasses)
                                        # read in each successive line from the interim release


# Flip SNPs to affy's order
m = match(newSnps2,ps2snpAll$AffySNPID)
toFlip = ps2snpAll$Reference[m]!=ps2snpAll$AlleleB[m]
flipThese = newData[,2]%in%newSnps[toFlip]
newData2 = newData
myFlipMat = newData2[flipThese,-c(1,2)]
myFlipMat[myFlipMat==0] = 5
myFlipMat[myFlipMat==2] = 0
myFlipMat[myFlipMat==5] = 2
newData2[flipThese,-c(1,2)] = myFlipMat

newDataOrdered = newData2[,-c(1:6)][snpOrder,]
newDataOrderedSNPs = newData2$SNP[snpOrder]

print( paste0(sum(flipThese) ," snps flipped back to affy allele order.") )


###### classify snps by MAF
print("Classifying SNPs by MAF (in all data, final release). ")

freqs = read.genotyped.maf(mafBins = c(0,1/1000,1/100,5/100,Inf))
                                        # flip minor alleles to match interim release
freqs$MinorAllele[ (freqs$SNP%in%newSnps[toFlip])&(freqs$freq>0.5) ] = 1
freqs$MinorAllele[ (freqs$SNP%in%newSnps[toFlip])&(!freqs$freq>0.5) ] = 2

OneisMinorAllele = newDataOrderedSNPs%in%freqs$SNP[freqs$MinorAllele==1]



# TO HERE

###### Compute differences

print("Counting concordances across interim samples in batch...")

#close(con)
con <- file(args[3], "r")
oldData = readLines(con,n=1)

cs = sapply(1:length(interimFam),function(i){

   # cs = sapply(1:3,function(i){
#for(i in 1:20){
    if(i%%100==0) print(i)
    iNew = which(newFam==interimFam[i])
    oldData = readLines(con,n=1)
    if( length(iNew)==0 ) return(rep(NA,5+5*2*(length(mafBins)-1)))
    dOld = str_split(oldData," ")[[1]][-c(1,2)]
    dOld = as.numeric(dOld[snpOrderOld])
    dOld[dOld=="-1"]=NA
    dNew = newDataOrdered[,iNew]
    
    c1 = sum(dOld!=dNew,na.rm=TRUE) # sum discordance excluding any missing
    c2 = sum((is.na(dOld))&(!is.na(dNew))) # sum changed from missing to non-missing
    c3 = sum((!is.na(dOld))&(is.na(dNew))) # sum changed from non-missing to missing
    nm = sum((!is.na(dOld))&(!is.na(dNew)))
    bm = sum((is.na(dOld))&(is.na(dNew)))
    all = c(c1,c2,c3,nm,bm)
                                        #   print(c1)
                                        # by frequency bin
    # snps where new data carries at least one copy of the minor allele.
    isMinor = ( ( dNew%in%c(1,2) )&( OneisMinorAllele ) )|( ( dNew%in%c(0,1) )&( !OneisMinorAllele ) )
    
    byBins = all
    byBins2 = c()
    for(b in 1:(length(mafBins)-1)){
        
        snps=freqs$SNP[freqs$bin==b]
        inBin = (newDataOrderedSNPs%in%snps)  # must carry minor allele in new data.
        inBinAndMinor = inBin&isMinor
        
        c1b = sum((dOld!=dNew)&(inBin),na.rm=TRUE) # sum discordance excluding any missing
        c2b = sum((is.na(dOld))&(!is.na(dNew))&(inBin)) # sum changed from missing to non-missing
        c3b = sum((!is.na(dOld))&(is.na(dNew))&(inBin)) # sum changed from non-missing to missing
        nmb = sum((!is.na(dOld))&(!is.na(dNew))&(inBin))
        bmb = sum((is.na(dOld))&(is.na(dNew))&(inBin))
        
        byBins =  c(byBins,c1b,c2b,c3b,nmb,bmb)
        
        c1b = sum((dOld!=dNew)&(inBinAndMinor),na.rm=TRUE) # sum discordance excluding any missing
        c2b = sum((is.na(dOld))&(!is.na(dNew))&(inBinAndMinor)) # sum changed from missing to non-missing
        c3b = sum((!is.na(dOld))&(is.na(dNew))&(inBinAndMinor)) # sum changed from non-missing to missing
        nmb = sum((!is.na(dOld))&(!is.na(dNew))&(inBinAndMinor))
        bmb = sum((is.na(dOld))&(is.na(dNew))&(inBinAndMinor))

        byBins2 =  c(byBins2,c1b,c2b,c3b,nmb,bmb)

    }
    
    return(c(byBins,byBins2))
})

close(con)


cs = t(cs)
basicColNames=c("nCallChanges","mToNm","NmTom","nonMissingBoth","MissingInBoth")
colnames(cs) = c(as.vector(t(sapply(basicColNames,function(x) paste0(x,c("",mafBins[-length(mafBins)]))))),
            as.vector(t(sapply(basicColNames,function(x) paste0(x,".Minor.",c(mafBins[-length(mafBins)]))))) )

rownames(cs) = interimFam

save(cs,intersectSnps,intersectSamples,file=outFile )


print( paste0("DONE! RData file saved: ",outFile) )
