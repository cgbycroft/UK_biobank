
h = c("/well/ukbiobank/expt/V2_QCed.SNP-QC/src/V2_QCed.snpqc-tests.R","/well/ukbiobank/expt/V2_QCed.SNP-QC/src/V2_QCed.bin2clusterplots.R","/well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/R/scripts/readPSperformance.R","/well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/R/scripts/auxFunctions.R")
for(s in h) source(s)




######## Which chromosome and what data are we plotting?

args = commandArgs(TRUE)

#args = c("test.out","1","plots", "-ymax","50","-hits", "/well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/GWAS/otherGWAS/GWAScatalogue/hg19/gwasCatalog-subset-Standing.height.RData", "-title", "-Euro-hits")

#args = c("/well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/otherGWAS/Intra.ocular.pressure..Goldmann.correlated..mean/BOLTLMM.v16/Intra.ocular.pressure..Goldmann.correlated..mean-BOLT-LMM-v16.out","genome","plots", "-ymax","50","-qc","-hits", "/well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/GWAS/otherGWAS/GWAScatalogue/hg19/gwasCatalog-subset-Intra.ocular.pressure..Goldmann.correlated..mean.RData", "-title", "-raw-hits")

print(args)

dataFile = args[1]
chroms = args[2]
plotOutDir = args[3]
if("-title"%in%args) extraTitle = args[which(args=="-title")+1] else extraTitle=""
useLmmInf=TRUE
if("-lreg"%in%args){
    useLmmInf=FALSE # use linear regression stats instead of lmm results
    extraTitle = paste0(extraTitle,"-lreg")
}

######## Get GWAS catalogue information
# NOTE: field descriptons are here: http://genome.ucsc.edu/cgi-bin/hgTables
columns = c("bin","chrom",
"chromStart",
"chromEnd",
"name",
"pubMedID",
"author",
"pubDate",
"journal",
"title",
"trait",
"initSample",
"replSample",
"region",
"genes",
"riskAllele",
"riskAlFreq",
"pValue",
"pValueDesc",
"orOrBeta",
"ci95",
"platform",
"cnv")

catFile = NULL
if("-hits"%in%args) {
    # NOTE: this overrides the -hi for QC colours
    catInfo = args[which(args=="-hits")+1]
    load(catInfo,verbose=TRUE)
    
    catFile = catPheno  # any ancestry
    colnames(catFile)[ncol(catFile)] = "Ancestry"
        
    catFile$Pvalue = catFile$V18
    catFile$SNP = catFile$V5
    catFile$BP = catFile$V4 # this is the chromEnd field
    catFile$CHR = gsub("chr","",catFile$V2)
    catFile$CHR[catFile$CHR%in%names(sexChroms)] = sexChroms[catFile$CHR[catFile$CHR%in%names(sexChroms)]]
    catFile$CHR = as.numeric(catFile$CHR)

    # adjust any Xchrom hits for PAR regions
    parHits = ( catFile$CHR == 23 ) & ( ( catFile$BP < par1 ) | ( catFile$BP > par2 ) )
    print( paste0( sum(parHits)," hits in PAR X") )
    catFile$CHR[parHits] = 25
    
    print( head(catFile) )
}


                                        # do we fix the y-axis?
Ymax = NULL
if("-ymax"%in%args) Ymax = as.numeric(args[which(args=="-ymax")+1])

                                        # which chroms?
if(chroms!="genome") {
    if(chroms=="all") {
        
        chroms = 1:22
        if( "-sex" %in% args ) chroms = 23:26
        
    } else {
        
        chroms = parse.range.string(chroms)        
    }
}


                                        # do we apply qc thresholds?
if("-qc" %in% args){
                                        # QC THRESHOLDS 
    minmaf = 0.001
    mininfo = 0.3
    maxmiss = 0.05  # maximum 5% missing data

} else {
    minmaf=0
    mininfo=0
    maxmiss=1
}


if("-minmaf" %in% args) minmaf = as.numeric(args[which(args=="-minmaf")+1])
if("-mininfo" %in% args) mininfo = as.numeric(args[which(args=="-mininfo")+1])
if("-maxmiss" %in% args) maxmiss = as.numeric(args[which(args=="-maxmiss")+1])


print("You asked for chromosomes")
print( chroms )

if(!grepl("%%",dataFile)){
                                        # read in the whole set (could be multiple chromosomes)
    Data = read.gwas(dataFile,chrom="genome",minmaf,mininfo,maxmiss,Ymax,QCexclusions=c(),extraTitle,useLmmInf=useLmmInf)
    rawPvalset=Data$Pvalset
                                        # are all the chromosomes in this dataset?
    chromsData = unique(Data$DF$CHR)
    if(chroms[1]!="genome") chroms = chroms[chroms%in%chromsData]
}

# colours for catalogue hits
traitCols = primaryColors[!primaryColors%in%c("blue","green")][1:length(table(catFile$V11))]
names(traitCols)=names(table(catFile$V11))

print("printing the following chromosomes")
print( chroms )



###################################
# Read in set of snps used in phasing!
phasedBims = read.phased.snps(chroms)

#### colours for manhattans
phasedColour = "red"
rareColour = "blue"
rareThreshold = 0.001


for(chrom in chroms){
    print(chrom)
    phased = phasedBims[[as.character(chrom)]]
                                        #phased$SNPimp = paste0(phased$V1,":",phased$V4,"_",phased$V5,"_",phased$V6)  # reconstruct ID in imputed data. Assume that alleles are in the same order!
    
    print(chrom)

    if(grepl("%%",dataFile)) {
        Data = read.gwas(gsub("%%",chrom,dataFile),chrom,minmaf,mininfo,maxmiss,Ymax,QCexclusions=c(),extraTitle,useLmmInf=useLmmInf)
    } 

    if(chrom=="genome") cexm=1.5 else cexm=1
    
    plot.BOLT.pvalues(Data,chrom,plotOutDir,plotQQ=FALSE,catFile=NULL,Ymax,moreToPlot=TRUE,extraFilename=paste0(".colourPhasedandRareLT",rareThreshold),cexManhattan=cexm)

    if(chrom!="genome"){
        
        inPhased = Data$DF$SNP%in%phased$V2
        points(Data$DF$BP[inPhased],-log10(Data$DF$P2[inPhased]),col=phasedColour,pch=16,cex=cexm)
                                        # rare snps (use colour of points border)
        rare = (Data$DF$MAF < rareThreshold) & (Data$DF$CHR==chrom)
        points(Data$DF$BP[rare],-log10(Data$DF$P2[rare]),col=rareColour,pch=1,lwd=2,cex=cexm)

        legend("topright",legend=c("Not phased","Phased",paste0("MAF < ",rareThreshold)),col=c("black",phasedColour,rareColour),bty="n",pch=c(16,16,1),cex=3,pt.lwd=c(1,1,2))
    }
                                        # add in the catalogue hits.
    if(!is.null(catFile)){
        print( "Printing catalogue hits..." )

        if(chrom=="genome"){
            
            catFile$BPorig=catFile$BP
            lastbase = 0
            
            for(i in unique(Data$DF$CHR) ){
                catFile$BP[catFile$CHR==i] = lastbase + catFile$BPorig[catFile$CHR==i]
                lastbase = lastbase + max(Data$DF$BP[Data$DF$CHR==i])
            }
            print(paste0(sum(is.na(catFile$BP))," positions in catalogue with no plotting index."))
            
        }

        if(chrom!="genome") catFileSub = catFile[catFile$CHR %in% chrom,] else catFileSub=catFile

        print( paste0( sum(-log10(catFileSub$Pvalue)>Ymax)," catalogue hits above ",Ymax) )
        
        catFileSub$Pvalue[-log10(catFileSub$Pvalue)>Ymax] = 10^(-Ymax)
                                        # plot non-European hits differently
        colors = traitCols[catFileSub$V11]
                                        #        colors[!grepl("European",catFileSub$Ancestry)] = "blue"
        print( table(catFileSub$Ancestry[!grepl("European",catFileSub$Ancestry)]) )
        
                                        # do we have these SNPs in UKBiobank? match on chrom and position
                                        #inHere = (catFileSub$BP %in% DF$BP)&(catFileSub$CHR == DF$CHR)
                                        #catFileSub$Pvalue[inHere] = DF$Pvalue[]
        
        points(catFileSub$BP,-log10(catFileSub$Pvalue),pch=8,col=colors,cex=4,lwd=1,xpd=NA)
        points(catFileSub$BP,-log10(catFileSub$Pvalue),pch=16,col=colors,cex=2,xpd=NA)
    }
    
    dev.off()
    
}
                                        # plot a separate legend plot
png(paste(plotOutDir,"/",gsub("chr*\\.m","chrgenome.m",Data$Pvalset),extraTitle,"-manhattan-LEGEND.png",sep=""),width=12,height=12,units="in",res=150,bg="transparent")
plot(NULL,xlim=c(0,1),ylim=c(0,1),axes=FALSE,xlab=NA,ylab=NA)
legend("top",horiz=FALSE,legend=names(traitCols),col=traitCols,pch=8,pt.cex=4,bty="n",cex=2,y.intersp=1.2)
legend("top",horiz=FALSE,legend=names(traitCols),col=traitCols,pch=16,pt.cex=2,bty="n",cex=2,text.col="transparent",y.intersp=1.2)
dev.off()

print(paste0("SEHR GUT! Plots saved in files with this kind of name: ",Data$Pvalset))


