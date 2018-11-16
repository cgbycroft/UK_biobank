

####  Extract intensities for interim release

#h = read.table('../qcoutput/Axiom_UKBL_Chip/ps2snp.txt',header=TRUE,stringsAsFactors=F)
h = read.table('../../qcoutput/Axiom_UKBL_Chip/sexchrom.txt',header=TRUE,stringsAsFactors=F)
h2 = read.table('../../qcoutput/Axiom_UKBB_Chip/sexchrom.txt',header=TRUE,stringsAsFactors=F)


source('/well/ukbiobank/qcoutput/QC-Scripts/R/scripts/bin2clusterplots.R')
source("/well/ukbiobank/qcoutput/QC-Scripts/R/scripts/batch2sampleinfo.R")

snp="Affx-34461684"
#h[h$AffySNPID==snp,]
snp = "Affx-79381704"  # <=== A UK Biobank only snp. Has two probesets! Best probset in each batch is: AX-94364249 or AX-94380998. In Interim release it was mostly AX-94364249. But in the final release it was: AX-94380998

#h2[h2$AffySNPID==snp,]
# 310 unique MT snps in interim; 253 have two probesets

pid = which(h2$AffySNPID==snp)[1]
pname = h2$ProbeSetID[pid]
    
#Batch = "UKBiLEVEAX_b11"
#Batch="UKBiLEVEAX_b5"
Batch="Batch_b010"
datadir = paste0("/well/ukbiobank/qcoutput/",Batch,"/GT1")
phenodir = paste0("/well/ukbiobank/qcoutput/",Batch)
phenofile = paste0(Batch,"_Sample_Table_Pheno.csv")
CsvFile = paste(phenodir,"/",phenofile,sep="")
CsvFile = paste0("/well/ukbiobank/qcoutput/",Batch,"/UKB_WCSGAX_",get.batch.no(Batch),"_Sample_Table_Pheno.csv")


data1 = unpack.probeset.data(datadir,pid=pid,is.autosomal=FALSE,pname=pname)
BatchInfo = get.batch.info(Batch,CsvFile)


png(file=paste0("clusterplot-example-MT-",snp,".png"),
    width=7,height=7,units="in",res=150)
cluster.plot(data1,main.text=Batch)
dev.off( )



####  Extract intensities for final release
library(stringr)
h = c("/well/ukbiobank/expt/V2_QCed.SNP-QC/src/V2_QCed.snpqc-tests.R","/well/ukbiobank/expt/V2_QCed.SNP-QC/src/V2_QCed.bin2clusterplots.R","/well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/R/scripts/readPSperformance.R","/well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/R/scripts/auxFunctions.R")
for(s in h) source(s)
source("/well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/R/scripts/ukbbcolors.R")


source('/well/ukbiobank/expt/V2_QCed.SNP-QC/src/V2_QCed.batch2sampleinfo.R')
source('/well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/R/scripts/bin2clusterplots_v2.R')
                                        # NOTE: location of calls (callsFiles) is defined in bin2clusterplots_v2.R
getSnpInfo <- function(SNPs,ps2snpTableAutosome,ps2snpTableSexchrom){

                                        # this reads the *.bim files for the order of the intensity/calls files.
    snpInfoAutosome = system( paste0("grep -E -wn \'",paste(SNPs,collapse="|"),"\' ",ps2snpTableAutosome) , intern=T)
    snpInfoSex = system( paste0("grep -E -wn \'",paste(SNPs,collapse="|"),"\' ",ps2snpTableSexchrom) , intern=T)

    snpInfo = c(snpInfoAutosome,snpInfoSex)
    
    allSnpInfo = t(sapply(snpInfo, function(x){    
        n = str_split(x,":")[[1]]
        p = str_split(n[2]," ")[[1]]
        c(as.numeric(n[1]) - 1,p)
    },simplify=T))    
    
    headerStuff = c("pid",read.table(ps2snpTableSexchrom,nrow=1,stringsAsFactors=F))
    
    if( length(allSnpInfo) == 0){
        allSnpInfo = as.data.frame(matrix(vector(),0,length(headerStuff) ))        
    } else {
        rownames(allSnpInfo) = 1:nrow(allSnpInfo)
    }
    colnames(allSnpInfo) = headerStuff

                                        # get pid from *bim files.
                                        # check that all SNPs are represented here
    if(sum(!SNPs%in%c(allSnpInfo[,"AffySNPID"],allSnpInfo[,"dbSNPRSID"])) > 0 ) {
        
        ids = SNPs[!SNPs%in%c(allSnpInfo[,"AffySNPID"],allSnpInfo[,"dbSNPRSID"])]
        
        print(paste0("WARNING: Unable to find markers ",paste(ids,collapse=", ")," in the axiom chip tables ",ps2snpTableAutosome," or ",ps2snpTableSexchrom,". Check that you have specified the right chip file. Will obviously not plot these markers (in some batches)."))
        
    }

    bimInfo = getPID(SNPs=allSnpInfo[,"dbSNPRSID"]) # function is defined in bin2clusterplots_v2.R
    print(head(bimInfo))
    print(head(allSnpInfo))

    allSnpInfo = cbind(allSnpInfo,bimInfo[match(allSnpInfo[,"dbSNPRSID"],bimInfo[,"pname"]),,drop=FALSE])
    
    return(allSnpInfo)
}


SNPs = snp
inputdir= "/well/ukbiobank/expt/V2_QCed.SNP-QC/data"
ps2snpTableAutosome = paste0(inputdir,"/V2_QCed.Axiom_UKBB_Chip/autosome.txt")
ps2snpTableSexchrom = paste0(inputdir,"/V2_QCed.Axiom_UKBB_Chip/sexchrom.txt")

UKBB = getSnpInfo(SNPs,ps2snpTableAutosome,ps2snpTableSexchrom)

ps2snpTableAutosome = paste0(inputdir,"/V2_QCed.Axiom_UKBL_Chip/autosome.txt")
ps2snpTableSexchrom = paste0(inputdir,"/V2_QCed.Axiom_UKBL_Chip/sexchrom.txt")

UKBL = getSnpInfo(SNPs,ps2snpTableAutosome,ps2snpTableSexchrom)

allSnpInfo = merge(UKBB,UKBL,by=colnames(UKBL)[-1],suffixes=c(".UKBB",".UKBL"),all=T,stringsAsFactors=F)

allSnpInfo = apply(allSnpInfo,2,as.character)

if(class(allSnpInfo)=="character") allSnpInfo = as.data.frame(t(allSnpInfo),stringsAsFactors=FALSE)

print(head(allSnpInfo))

if( sum(!is.na(allSnpInfo[,"line"]))==0 ) {
    print("WARNING: No data found for any snps in release! Stopping.")
    quit()
}

print("Reading in sample qc information based on released data...")
releaseSampleQCPrefix = 'b1__b11-b001__b095-sampleTable_v4'
load(paste0('/well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/ForRelease/',releaseSampleQCPrefix,'_allColumns.RData'),verbose=TRUE)
                                        #    batchInfo = read.delim(paste0('/well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/ForRelease/',releaseSampleQCPrefix,'_allColumns.csv'),header=TRUE,stringsAsFactors=FALSE,sep=",")

                                        # Read in order of samples in fam files.
famOrder = read.table(paste0(callsFiles,".chr22.fam"),header=FALSE,stringsAsFactors=FALSE)[,2]
outTable = outTable[match(famOrder,outTable$PIID), ]

##########
i=1
##########

snpName = allSnpInfo[i,"AffySNPID"]
pname = allSnpInfo[i,"ProbeSetID"]
rsID = allSnpInfo[i,"dbSNPRSID"]
chrom = allSnpInfo[i,"Chromosome"]
position = as.numeric(allSnpInfo[i,"Position"])
isPAR = sum(as.numeric(allSnpInfo[i,c("IsInPAR1","IsInPAR2")]))
bim = allSnpInfo[i,"bim"]

pid.bim = as.numeric(allSnpInfo[i,"line"])

bimChrom = str_split(str_split(basename(bim),"chr")[[1]][2],"\\.")[[1]][1]

calls = unpack.calls2(paste0(callsFiles,".chr",bimChrom,".bed"),famFile=NULL,pid.bim)
intensities = unpack.intensities2(paste0(intensityFiles,".chr",bimChrom,".bin"),famFile=NULL,pid.bim)

dataAll = list(calls=calls,intensities=intensities)

print( paste0(callsFiles,".chr",bimChrom,".bed") )
print( paste0(intensityFiles,".chr",bimChrom,".bin") )
print(rsID)

subset = (outTable$Batch==Batch)&(!is.na(outTable$Batch))
batchInfo2 = outTable[subset,] # only plot samples in this batch. outTable should be in the order of famOrder.

n = length(intensities)/2
logratio2All = intensities[2*(1:n)-1]
strength2All = intensities[2*(1:n)]
names(logratio2All) = names(strength2All) = outTable$Best.Array
a2all = 2^(strength2All + logratio2All/2)
b2all = 2^(strength2All - logratio2All/2)

### By batch
intensitySubset = sort(c(which(subset)*2,which(subset)*2-1))
data = list(calls=dataAll$calls[subset],intensities=dataAll$intensities[intensitySubset])

data2 = data


Rscript /well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/R/scripts/V2_QCed.makeClusterPlots_v2.R /well/ukbiobank/expt/V2_QCed.SNP-QC/data Affx-34461684 cluster_plots_MT -b UKBiLEVEAX_b1,UKBiLEVEAX_b3,UKBiLEVEAX_b5,Batch_b001,Batch_b008,Batch_b016 -rNumber 8194 &

Rscript /well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/R/scripts/V2_QCed.makeClusterPlots_v2.R /well/ukbiobank/expt/V2_QCed.SNP-QC/data Affx-79381704 cluster_plots_MT -b UKBiLEVEAX_b1,UKBiLEVEAX_b3,UKBiLEVEAX_b5,Batch_b001,Batch_b008,Batch_b016 -rNumber 8195 &

    
################ Compare
IND="A550465-4276624-031217-258_C03"
sampleName=outTable$Sample.Name[outTable$PIID==IND]

# Was in UKBileve 11 in final data
                                        # Not
for(Batch in all.batches()[1:33]){
    Batch1=Batch
    if(Batch%in%ukbiobank.batches()) Batch1 = paste0("UKB_WCSGAX_",get.batch.no(Batch))
    CsvFile = paste0("/well/ukbiobank/qcoutput/",Batch,"/",Batch1,"_Sample_Table_Pheno.csv")
    b = get.batch.info(Batch,CsvFile)
    if(sampleName%in%b$Sample.Name) print(Batch)
}


int1 = data1$intensities
n = length(int1)/2
logratio1 = int1[2*(1:n)-1]
strength1 = int1[2*(1:n)]
names(logratio1) = names(strength1) = BatchInfo$Best.Array
a1 = 2^(strength1 + logratio1/2)
b1 = 2^(strength1 - logratio1/2)


int2 = data2$intensities
n = length(int2)/2
logratio2 = int2[2*(1:n)-1]
strength2 = int2[2*(1:n)]
names(logratio2) = names(strength2) = batchInfo2$PIID
a2 = 2^(strength2 + logratio2/2)
b2 = 2^(strength2 - logratio2/2)


cal1 = data1$calls
names(cal1) = BatchInfo$Best.Array

cal2 = data2$calls
names(cal2) = batchInfo2$PIID

hist(logratio1[names(logratio2)]/logratio2,plot=F)

hist(cal1[names(cal2)]/cal2,plot=F)

diffs = abs(a1[names(a2)]-a2)
hist(diffs,plot=F)

png(paste0("intensities-",snp,"_A.png"),height=1000,width=1000,res=150)
plot(log(a1[names(a2)]),log(a2),main=Batch)
abline(0,1,col="red")
dev.off()

png(paste0("intensities-",snp,"_B.png"),height=1000,width=1000,res=150)
plot(log(b1[names(b2)]),log(b2),main=Batch,
     xlab="log( B intensity ) interim data",
     ylab="log( B intensity ) final data")
abline(0,1,col="red")
dev.off()




####### SNP info in final release

BB.ps2snp = ukbiobank.ps2snp("sexchrom")
BL.ps2snp = ukbileve.ps2snp("sexchrom")
ps2snpBoth = unique(BL.ps2snp[BL.ps2snp$AffySNPID%in%BB.ps2snp$AffySNPID,])
ps2snpEither = unique(rbind(BL.ps2snp,BB.ps2snp))

sum(ps2snpEither$Chromosome==26)
# 266 MT snps in final release

sum(ps2snpEither$AffySNPID[ps2snpEither$Chromosome==26]%in%both$AffySNPID)
# all are in previous release

mtsnps = table(both$AffySNPID[both$Chromosome=="MT"])
twoPobes = names(mtsnps)[mtsnps>1]
sum(ps2snpEither$AffySNPID[ps2snpEither$Chromosome==26]%in%twoPobes)

# 214 have two probesets
