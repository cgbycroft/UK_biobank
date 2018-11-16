#########
# plot ychrom depletion stats
#########

setwd('/well/donnelly/ukbiobank_project_8874/clare/Ychrom')

#############
vers="v1"   #
#############


# All samples
genoInfoFile = '/well/donnelly/ukbiobank_project_8874/clare/phenotypes/ukb8874-basic-allGenotypedSamples-processed.txt'  # this file was created by running clare/commonScripts/getCSVinfo-ukb8874.py, and it has all the phenotypes I want.

genoInfo = read.delim(genoInfoFile,header=TRUE,stringsAsFactors=FALSE,skip=1,sep=",")


# Samples used in this version of the GWAS. Created in subsetPhenotypes.R
genoInfoGWAS = read.table(paste0("data/YchromPhenotypesForBOLT-",vers,".txt"),stringsAsFactors=FALSE,header=TRUE)


var1 = "Average.Y.chromosome.intensities.for.determining.sex"
var2 = "Age.when.attended.assessment.centre"
var3 = "Average.X.chromosome.intensities.for.determining.sex"

males = (genoInfo$Genetic.sex==1)&(genoInfo$Sex==1)
keep = is.na(genoInfo$Recommended.genomic.analysis.exclusions)
sampled = genoInfo$Encoded.anonymised.participant.ID%in%genoInfoGWAS$IID

x = genoInfo[males&keep,var1]
png(paste0("scripts/depletion/plots/YchromIntensity-hist-allmales.png"),height=1000,width=1000,res=150)
hist(x,breaks=1000,xlab="Mean Y chromosome intensity")
dev.off()

x = genoInfo[males&keep,var1]/genoInfo[males&keep,var3]
png(paste0("scripts/depletion/plots/XYIntensityRatio-hist-allmales.png"),height=1000,width=1000,res=150)
hist(x,breaks=1000,xlab="Y:X intensity ratio")
dev.off()

x = genoInfo[males&keep&sampled,var1]
png(paste0("scripts/depletion/plots/YchromIntensity-hist-WhiteBritishMales.png"),height=1000,width=1000,res=150)
hist(x,breaks=1000,xlab="Mean Y chromosome intensity")
dev.off()

x = genoInfo[(!males)&keep,var1]
png(paste0("scripts/depletion/plots/YchromIntensity-hist-allfemales.png"),height=1000,width=1000,res=150)
hist(x,breaks=1000,xlab="Mean Y chromosome intensity")
dev.off()


x = genoInfo[keep,var1]
y = genoInfo[keep,var3]
colors=rep("red",nrow(genoInfo))
colors[males]="blue"

png(paste0("scripts/depletion/plots/YchromIntensity-byXchromIntensity-allSamples.png"),height=1000,width=1000,res=150)
plot(x,y,col=colors[keep],xlab="Y intensity",ylab="X intensity")
dev.off()



# 2-year age-groups
png(paste0("scripts/depletion/plots/YchromIntensity-hist-allmales-byAge-%02d.png"),height=1000,width=1000,res=150)
par(mfrow=c(3,1),mar=c(2,3,1,1))
x = genoInfo[males&keep,var1]
xlims=c(min(x),max(x))
for(a in seq(40,69,2)){
    age = genoInfo$Age.when.attended.assessment.centre%in%c(a,a+1)
    x = genoInfo[males&keep&age,var1]
    hist(x,breaks=500,xlim=xlims,xlab="Mean Y chromosome intensity",main="",border=NA,col="black")
    legend("topright",bty="n",legend=paste(a,a+1),cex=2,inset=c(-2,0))
}

dev.off()


png(paste0("scripts/depletion/plots/YchromIntensity-hist-allmales-byAge-zoom-%02d.png"),height=1000,width=1000,res=150)
par(mfrow=c(3,1),mar=c(2,3,1,1))
x = genoInfo[males&keep,var1]
xlims=c(min(x),500)
for(a in seq(40,69,2)){
    age = genoInfo$Age.when.attended.assessment.centre%in%c(a,a+1)
    x = genoInfo[males&keep&age,var1]
    hist(x,breaks=500,xlim=xlims,xlab="Mean Y chromosome intensity",main="",border=NA,col="black")
    legend("topright",bty="n",legend=paste(a,a+1),cex=2,inset=c(-2,0))
}

dev.off()


png(paste0("scripts/depletion/plots/XYIntensityRatio-hist-allmales-byAge-zoom-%02d.png"),height=1000,width=1000,res=150)
par(mfrow=c(3,1),mar=c(2,3,1,1))
x = genoInfo[males&keep,var1]/genoInfo[males&keep,var3]
xlims=c(min(x),max(x))
for(a in seq(40,69,2)){
    age = genoInfo$Age.when.attended.assessment.centre%in%c(a,a+1)
    x = genoInfo[males&keep&age,var1]/genoInfo[males&keep&age,var3]
    hist(x,breaks=500,xlim=xlims,xlab="Mean Y chromosome intensity",main="",border=NA,col="black")
    legend("topright",bty="n",legend=paste(a,a+1),cex=2,inset=c(-2,0))
}

dev.off()
