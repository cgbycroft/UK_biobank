library(zoo)
plink='/well/ukbiobank/qcoutput/Src/plink-v1.90b3k-x86_64/plink'

args=commandArgs(trailingOnly=T)

args = c('/well/ukbiobank/expt/V2_QCed.SNP-QC/src/V2_QCed.snpqc-tests.R','/well/ukbiobank/expt/V2_QCed.SNP-QC/src/V2_QCed.bin2clusterplots.R','/well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/R/scripts/readPSperformance.R','/well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/R/scripts/auxFunctions.R','-in','/well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/PCA/b1__b11-b001__b095-autosome-sampleqc-fastpca-highquality-round2-coolibah.snpload.map','-out','b1__b11-b001__b095-pca-UKbio-round2')


print(args)
h = args[-c(which(args%in%c("-in","-out")),1+which(args%in%c("-in","-out")))]
for(helperScript in h){
    source(helperScript)
}

OutputFile = args[which(args%in%c("-out"))+1]
OutDir =  paste0(baseSampleQCDir,"/QC-Scripts/PCA/pca-UKBio")
SnpLoads = args[which(args%in%c("-in"))+1]


# read in PCA data (created in plot-pca-UKBio.R)
Rfile=paste0(baseSampleQCDir,"/data/PCA/",OutputFile,".RData")
load(Rfile,verbose=TRUE)

nPCs = sum(grepl("PC",colnames(PCs)))

# read in snploads
loads = dplyr::tbl_df(read.table(SnpLoads,header=FALSE,stringsAsFactors=FALSE))

# read in other sample info
Info=read.multiple.batch.info(batchInfoFields)
Info = dplyr::tbl_df(Info)

# read in het-missing info
load(paste0(baseSampleQCDir,"/data/Combined/b1__b11-b001__b095-autosome-sampleqc-hetcorrected-6PCs-imiss.RData"),verbose=TRUE)
hetmissOutliers = read.table(paste0(baseSampleQCDir,"/QC-Scripts/Sample-QC/HetMissing/b1__b11-b001__b095-autosome-sampleqc-hetcorrected-6PCs-imiss-all-hetmiss-outliers.txt"),header=FALSE)[,1]


# read in kinship data
kin = read.table(paste0(baseSampleQCDir,"/data/Relatedness/b1__b11-b001__b095-pair_batches.filtered-hetmiss-out.kin0"),header=TRUE,stringsAsFactors=FALSE)
kin$class = get.kin.classes(kin,0.0012)

kinOrig = read.table(paste0(baseSampleQCDir,"/data/Relatedness/b1__b11-b001__b095-pair_batches.filtered.kin0"),header=TRUE,stringsAsFactors=FALSE)
kinOrig$class = get.kin.classes(kinOrig,0.0012)

# who was in the PC-calculation?
fam = read.table(paste0(baseSampleQCDir,"/data/PCA/b1__b11-b001__b095-autosome-sampleqc-fastpca-round2-highquality.fam"),header=FALSE,stringsAsFactors=FALSE)

PCs$inpcCalc = PCs$PIID%in%fam$V1
# 407599 samples

###########
# CHECK: PC 27-28 which has flying samples


# crudely select samples in top-right corner
pc27 = PCs[["PC27"]]
pc28 = PCs[["PC28"]]

is.in.ellipse <- function(x1,y1,x0,vx,y0,vy,cov,c=2){  
  theta = 0.5 * atan2(2*cov, vx-vy)
  sint = sin(theta)
  cost = cos(theta)
  a = c*sqrt(vx*cost*cost + vy*sint*sint + cov*2*sint*cost)
  b = c*sqrt(vx*sint*sint + vy*cost*cost - cov*2*sint*cost)
  sint = sin(theta)
  cost = cos(theta)
  test = ((cost*(x1 - x0) + sint*(y1 - y0))^2)/a^2 + ((sint*(x1 - x0) - cost*(y1 - y0))^2)/b^2
  return(test <= 1)
}

# For PC 27 and 28; define the ellipse
cov = 0 # the off-diagonals of the covariance matrix
x0=mean(pc27) # the centre of the ellipse
y0=mean(pc28)
vx=36 # the variances in x and y directions. The width of Ellipse will be 2*sqrt(vx)
vy=36

ellipsePoints = ellipse.info(x0,vx,y0,vy,cov)

#Order = sample(1:nrow(PCs),1000)
Order = order.by.number.occurrences(PCs[["Colors"]])

png(paste0(OutDir,'/plots/',OutputFile,'-pc27-28-ellipse.png'),width=1500,height=1500,res=150)
plot(PCs[["PC27"]][Order],PCs[["PC28"]][Order],col=PCs$Colors[Order],pch=PCs$Chars[Order],xlab="PC27",ylab="PC28")
lines(ellipsePoints[,1],ellipsePoints[,2],col="red",lty=3,lwd=2)
dev.off()

isIn = is.in.ellipse(PCs[["PC27"]],PCs[["PC28"]],x0,vx,y0,vy,cov)

toCheck = (!isIn)&(PCs[["PC27"]]>0)&(PCs[["PC28"]]>0)

cols = PCs$Colors
cols[!toCheck] = "gray"
png(paste0(OutDir,'/plots/',OutputFile,'-pc27-28-SamplesToCheck.png'),width=1500,height=1500,res=150)
plot(PCs[["PC27"]][Order],PCs[["PC28"]][Order],col=cols[Order],pch=PCs$Chars[Order],xlab="PC27",ylab="PC28")
lines(ellipsePoints[,1],ellipsePoints[,2],col="red",lty=3,lwd=2)
dev.off()


#####
samplesToCheck = PCs$PIID[toCheck]
#####

######
# Look at properties of these samples

# hetmissing
y = Table$het.corrected
x = Table$logit.miss
cols = rep("grey",length(x))
cols[Table$IID%in%samplesToCheck] = "red"
# How many were actually used in the PC calculation?
cols[(Table$IID%in%samplesToCheck)&(Table$IID%in%PCs$PIID[PCs$inpcCalc])] = "blue"

Order = order.by.number.occurrences(cols)
png(paste0(OutDir,'/plots/',OutputFile,'-pc27-28-SamplesToCheck-hetmiss.png'),width=1500,height=1500,res=150)
plot(x[Order],y[Order],col=cols[Order])
dev.off()


# Are they the samples that had the artefact????
artefact = read.table("/well/ukbiobank/qcoutput.V2_QCed.sample-QC-testing/QC-Scripts/PCA/pca-UKBio/b1__b11-b001__b095-pca-UKbio-proj40-dim50-pc17-19-samplesInCluster.txt",header=FALSE)[,1]

sum(samplesToCheck%in%artefact)
# there is only one sample!!

# Are they the samples that weird stuff in pc17 last time????
pc19samples = read.table("/well/ukbiobank/qcoutput.V2_QCed.sample-QC-testing/QC-Scripts/PCA/pca-UKBio/b1__b11-b001__b095-pca-UKbio-proj40-dim50-pc19-samplesWithScoreGt20.txt",header=FALSE)[,1]

sum(samplesToCheck%in%pc19samples)
# almost all of them!!

# How many 3rd-degree relatives do they have (after excluding hetmiss outliers)?
kinTotal = table(c(kin$ID1[kin$class=="3rd degree"],kin$ID2[kin$class=="3rd degree"]))
range(kinTotal[samplesToCheck],na.rm=TRUE)

kinTotalO = table(c(kinOrig$ID1[kinOrig$class=="3rd degree"],kinOrig$ID2[kinOrig$class=="3rd degree"]))
range(kinTotalO[samplesToCheck],na.rm=TRUE)
b = samplesToCheck[samplesToCheck%in%PCs$PIID[PCs$inpcCalc]]
range(kinTotalO[b],na.rm=TRUE)
# 1 to 38 relatives. 220 out of 251 samples have at least one relative.

# Any relationship with x/y intensity?? Contaminant samples might contain a bit of Y chromosome?
y = Info$Y.intensity/Info$X.intensity
x = pc27[match(Info$PIID,PCs$PIID)]
cols = rep("grey",length(x))
cols[Info$PIID%in%samplesToCheck] = "red"
# How many were actually used in the PC calculation?
cols[(Info$PIID%in%samplesToCheck)&(Info$PIID%in%fam$V2)] = "blue"
 
Order = order.by.number.occurrences(cols)
png(paste0(OutDir,'/plots/',OutputFile,'-pc27-28-SamplesToCheck-YXintensity.png'),width=1500,height=1500,res=150)
plot(x[Order],y[Order],col=cols[Order],xlab="PC27",ylab="Y:X intensity")
dev.off()

# Which batches are the samples in?
sort(table(Info$Batch[Info$PIID%in%samplesToCheck]))
sort(table(Info$Batch[Info$PIID%in%intersect(samplesToCheck,fam$V2)]))


# Which SNPs are strongly associated with these samples?
write.table(cbind(samplesToCheck,samplesToCheck),file=paste0(OutputFile,"-pc27-28-samplesInClusterPlink.txt"),col.names=FALSE,quote=FALSE,row.names=FALSE)
write.table(samplesToCheck,file=paste0(OutputFile,"-pc27-28-samplesInCluster.txt"),col.names=FALSE,quote=FALSE,row.names=FALSE)
write.table(intersect(samplesToCheck,fam$V2),file=paste0(OutputFile,"-pc27-28-samplesInClusterAndpcaCalculation.txt"),col.names=FALSE,quote=FALSE,row.names=FALSE)

system(paste0(plink," --bfile ../../../data/Combined/b1__b11-b001__b095-autosome-sampleqc --keep ",OutputFile,"-pc27-28-samplesInClusterPlink.txt --keep-allele-order --freqx --out ",OutputFile,"-pc27-28-cluster"))
system(paste0(plink," --bfile ../../../data/Combined/b1__b11-b001__b095-autosome-sampleqc --keep-allele-order --freqx --out b1__b11-b001__b095-autosome-sampleqc"))

counts1 = read.delim(paste0(OutputFile,"-pc27-28-cluster.frqx"),stringsAsFactors=F,header=T,sep="\t")
counts2 = read.delim(paste0("b1__b11-b001__b095-autosome-sampleqc.frqx"),stringsAsFactors=F,header=T,sep="\t")

counts1 = tbl_df(counts1)
counts2 = tbl_df(counts2)
counts = left_join(counts1,counts2,by=c("SNP"="SNP"))

gens1 = counts[,c(5:7,10)]
gens2 = counts[,c(14:16,19)]
toCompare=cbind(counts$CHR.x,(gens2-gens1)/2,(gens2-gens1)/2,gens1/2,gens1/2) # /2 just because routine assumes male/female are separate

# run fisher exact test on counts between group and non-group
compare = apply(toCompare,1,CompareGenotypes)
compare=compare[1,] # just the total comparison`
names(compare) = counts$SNP

# plot the worst offenders
offend = names(sort(compare))[1:50]

forSystem = paste0("Rscript /well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/R/scripts/V2_QCed.makeClusterPlots.R /well/ukbiobank/expt/V2_QCed.SNP-QC/data ",paste(offend,collapse=",")," /well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/PCA/pca-UKBio/plots/cluster_plots-pc27-28-samplesInClusterAndpcaCalculation-signifSNPs-orig -b Batch_b061,Batch_b064,Batch_b092,Batch_b095,UKBiLEVEAX_b7 -samples /well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/PCA/pca-UKBio/b1__b11-b001__b095-pca-UKbio-round2-pc27-28-samplesInClusterAndpcaCalculation.txt -orig")

system(forSystem)


# plot DQC vs sample CR
x = Info$dQC
y = Table$miss[match(Info$PIID,Table$IID)]
cols = rep("grey",length(x))
cols[Info$PIID%in%samplesToCheck] = "red"
# How many were actually used in the PC calculation?
cols[(Info$PIID%in%samplesToCheck)&(Info$PIID%in%fam$V2)] = "blue"
cols[(Info$PIID%in%hetmissOutliers)] = "green"
 
Order = order.by.number.occurrences(cols)

png(paste0(OutDir,'/plots/',OutputFile,'-pc27-28-SamplesToCheck-MissingbydQC.png'),width=1500,height=1500,res=150)
plot(x[Order],y[Order],col=cols[Order],xlab="dQC",ylab="Missing rate Sample QC")
dev.off()

y = Info$Cluster.CR
png(paste0(OutDir,'/plots/',OutputFile,'-pc27-28-SamplesToCheck-Cluster.CRbydQC.png'),width=1500,height=1500,res=150)
plot(x[Order],y[Order],col=cols[Order],xlab="dQC",ylab="Cluster.CR")
dev.off()
