## Script to plot (raw) heterozygosity and missingness
library(dplyr)

args = commandArgs(trailingOnly=T)
#args = c("/well/ukbiobank/expt/V2_QCed.SNP-QC/src/V2_QCed.snpqc-tests.R","/well/ukbiobank/expt/V2_QCed.SNP-QC/src/V2_QCed.bin2clusterplots.R","/well/ukbiobank/qcoutput.V2_QCed.sample-QC-testing/QC-Scripts/R/scripts/readPSperformance.R","/well/ukbiobank/qcoutput.V2_QCed.sample-QC-testing/QC-Scripts/R/scripts/auxFunctions.R",
#"-in","/well/ukbiobank/qcoutput.V2_QCed.sample-QC-testing/data/Combined/b1__b11-b001__b095-autosome-sampleqc","-npcs","6PCs")

print(args)

h = args[-c(which(args%in%c("-in","-out","-npcs")),1+which(args%in%c("-in","-out","-npcs")))]
for(helperScript in h){
    source(helperScript)
}

input = args[which(args=="-in") + 1]
nPCs = args[which(args=="-npcs") + 1]
outname = paste0("plots/",basename(input),"-",nPCs)  # just an adjustment for testing. %%%% change this back to the original simple name later.

# read in Rdata from correction of het
load(paste0(input,"-hetcorrected-",nPCs,"-imiss.RData"),verbose=T)

# read in LROH data (if you have it)

lrohFile = paste0(dirname(input),"/../Relatedness/",basename(input),"-1000.KB.hom.indiv")
if(basename(lrohFile) %in% list.files(dirname(lrohFile))){
    plotLROH=TRUE
    lroh = read.table(lrohFile,header=T,stringsAsFactors=F)
    lroh = dplyr::tbl_df(lroh)    
} else {
    plotLROH=FALSE
    print("NOTE: there is no lroh file... you better run ./compute-LROH.sh! Never mind, I'll carry on printing the heterozygosity stats. :D")
}
   
# get other info
moreInfo = read.multiple.batch.info(fields = c("Ethnic.background","Country.of.Birth..non.UK.origin.","Country.of.birth..UK.elsewhere."))
moreInfo$Country = get.place.of.birth(moreInfo)


# combine with Table
if(plotLROH) Table = left_join(Table,lroh,by=c("IID"="IID"))
Table = left_join(Table,moreInfo,by=c("IID"="PIID"))

#### exclude any reference samples
referenceList=paste0(baseSampleQCDir,"/QC-Scripts/referenceList.txt")
references = read.table(referenceList,header=F,stringsAsFactors=F)[,1]

nrefs = sum(Table$IID%in%references)

if(nrefs > 0 ){
    print(paste0("WARNING: ",nrefs," reference sample are found in het file ",hetFile,". Do you expect this?" ))
    print("Removing them from output in any case...")
    Table = Table[!(Table$IID%in%references),]
}


## get colours for country of birth, if PC-corrected het is lower than 0.18 (for lroh plots)
load(paste0(baseSampleQCDir,"/QC-Scripts/PCA/pca-UKBio/pca-UKBio-countryColours.RData"),verbose=T)

countries = sort(table(Table$Country[( Table$het.corrected < 0.18 ) ]),decreasing=T)
countriesNames = names(countries)[countries > 5]
h = generateColours(countriesNames)

Table$Colors2 = h[[1]][Table$Country]
Table$Chars2 = h[[2]][Table$Country]
Table$Chars2[is.na(Table$Chars2)] = 1
Table$Colors2[is.na(Table$Colors2)] = "black"


#### plot it

## ethnicity legend
eths = ethnicities[ethnicityOrder]
eths = eths[eths%in%unique(Table$Pops)]

png(paste0(outname,"-correctedHetbyMissing-Legend.png"),width=1000,height=1000,res=150)
plot.new()
legend("top",legend=eths,col=ethnicity2col[eths],pch=ethnicity2char[eths],bty="n",pt.lwd=2,horiz=F)
dev.off()


## Order the samples so that the most numerous ethnicity (British) is added to the scatter plot first
Order = order.by.number.occurrences(Table$Colors)

## plot raw heterozygosity by missing rate
y = Table$het[Order]
x = Table$logit.miss[Order]

png(paste0(outname,"-uncorrectedHetbyMissing.png"),width=1000,height=1000,res=150)
    plot(x,y,ylab="heterozygosity" ,xlab="logit missing",col=Table$Colors[Order],pch=Table$Chars[Order])
dev.off()

y = Table$het.corrected[Order]
png(paste0(outname,"-correctedHetbyMissing.png"),width=1000,height=1000,res=150)
    plot(x,y,ylab="PC-corrected heterozygosity" ,xlab="logit missing",col=Table$Colors[Order],pch=Table$Chars[Order])
dev.off()

x = Table$het[Order]
png(paste0(outname,"-correctedHetbyUncorrectedHet.png"),width=1000,height=1000,res=150)
    plot(x,y,xlab="Raw heterozygosity" ,ylab="PC-corrected heterozygosity",col=Table$Colors[Order],pch=Table$Chars[Order])
dev.off()


   
## plot LROH
if(plotLROH){
    y = Table$het.corrected[Order]
    x = Table$KB[Order]

    png(paste0(outname,"-correctedHetbyLROH-TotalKB.png"),width=1000,height=1000,res=150)
    plot(x,y,ylab="PC-corrected heterozygosity" ,xlab="Total length of ROH segments (KB)",col=Table$Colors[Order],pch=Table$Chars[Order])
    dev.off()

    png(paste0(outname,"-correctedHetbyLROH-TotalKB-COB.png"),width=1000,height=1000,res=150)
    plot(x,y,ylab="PC-corrected heterozygosity" ,xlab="Total length of ROH segments (KB)",col=Table$Colors2[Order],pch=Table$Chars2[Order])
    dev.off()


    png(paste0(outname,"-correctedHetbyLROH-TotalKB-COB-Legend.png"),width=1000,height=1000,res=150)
    plot.new()
    legend("topleft",legend=names(h[[1]])[1:20],col=h[[1]][1:20],pch=h[[2]][1:20],bty="n",pt.lwd=2,horiz=F)
    legend("top",legend=names(h[[1]])[21:length(h[[1]])],col=h[[1]][21:length(h[[1]])],pch=h[[2]][21:length(h[[1]])],bty="n",pt.lwd=2,horiz=F)
    dev.off()
    
}

## plot each ethnicity separately
png(paste0(outname,"-correctedHetbyMissing-byEthnicity.png"),width=2000,height=2000,res=150)
par(mfrow=c(4,4))    
ylims= c( min(Table$het.corrected), max(Table$het.corrected) )
xlims= c( min(Table$logit.miss), max(Table$logit.miss) )

for(eth in eths){
    
    print(eth)
    these = Table$Pops==eth
    x = Table$logit.miss[these]
    y = Table$het.corrected[these]
    
    plot(x,y,ylab="PC-corrected heterozygosity" ,xlab="logit missing",col=Table$Colors[these],pch=Table$Chars[these],
         main=paste0(eth," (",sum(these),")" ) ,ylim=ylims ,xlim=xlims)
    abline(mean(y),0,lty=3)
}

dev.off()


# for 'other/unknown' plot by country of birth
unknowns = Table[Table$Pops=="Other/Unknown",]
countries = sort(table(unknowns$Country),decreasing=T)


png(paste0(outname,"-correctedHetbyMissing-Unknown-byCOB-%02d.png"),width=2000,height=2000,res=150)
par(mfrow=c(4,4))    

for(country in names(countries)[countries > 30]){
    # countries with more than 20 samples
    print(country)
    these = unknowns$Country==country
    x = unknowns$logit.miss[these]
    y = unknowns$het.corrected[these]
    
plot(x,y,ylab="PC-corrected heterozygosity" ,xlab="logit missing",col=unknowns$Colors2[these],pch=unknowns$Chars2[these],
         main=paste0(country," (",sum(these),")" ) ,ylim=ylims ,xlim=xlims)
abline(mean(y),0,lty=3)
    
}

dev.off()

print("DONE!")
