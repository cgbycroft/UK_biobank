########
# usage:
# Rscript /well/ukbiobank/expt/V2_QCed.SNP-QC/src/V2_QCed.makeClusterPlots.R <input directory with batch-based binary files> <comma-separated list of Affy SNPIDs, or file name with list> <output directory> [-b <batchName, or range, or comma-separated list>] [-colour <variable in sample table. See colourOptions>] [-orig]
#
# DEFAULTS: prints snp at all batches, 9 batches to a page
########

library(stringr)

#source('/well/ukbiobank/expt/V2_QCed.SNP-QC/src/V2_QCed.bin2clusterplots_v2.R')
source('/well/ukbiobank/expt/V2_QCed.SNP-QC/src/V2_QCed.batch2sampleinfo.R')
source('/well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/R/scripts/bin2clusterplots_v2.R')

args = commandArgs(trailingOnly=T)
#args = c("/well/ukbiobank/expt/V2_QCed.SNP-QC/data","Affx-13890966,Affx-13890977,Affx-13890978,Affx-52187954,Affx-80256891","/well/ukbiobank/expt/V2_QCed.SNP-QC/data/plots/test_clusterPlots_plate","-b","Batch_b080-Batch_b085","-colour","Processed.Plate.Name","-samples","/well/ukbiobank/qcoutput.V2_QCed.sample-QC-testing/QC-Scripts/PCA/pca-UKBio/b1__b11-b001__b095-pca-UKbio-proj40-dim50-pc17-snpsWithHighloads0.05.txt")

# args = c("/well/ukbiobank/expt/V2_QCed.SNP-QC/data","rs3131962","test","-b","UKBiLEVEAX_b1,Batch_b052","-released")

print(args)

inputdir = args[1]
outputdir = args[3]

allBatches = all.batches()

if("-b"%in%args) {
    batches = args[which(args=="-b") + 1]
    batchesPlotted = paste(sapply(batches,get.batch.no),collapse="-")

    if(grepl("-",batches)) {
        batchRange = str_split(batches,"-")[[1]]
        batches = allBatches[which(allBatches==batchRange[1]):which(allBatches==batchRange[2])]
        batchesPlotted = paste(get.batch.no(batches[1]),get.batch.no(batches[length(batches)]),sep="--")    
    } else {
        if(grepl(",",batches)) {
            batches = str_split(batches,",")[[1]]
            batchesPlotted = paste(get.batch.no(batches[1]),get.batch.no(batches[length(batches)]),sep="--")
        }
    }
} else {
    batches = allBatches
    sapply(batches,get.batch.no)[1]
    batchesPlotted = paste(get.batch.no(batches[1]),get.batch.no(batches[length(batches)]),sep="--")    
}

batches = unique(batches)

SNPs = args[2]
if(grepl(",",SNPs)) SNPs = str_split(args[2],",")[[1]] else if(grepl(".txt",SNPs)) SNPs = read.table(args[2],stringsAsFactors=F,header=F)[,1]

# do we use the released data or the qc data?
if("-released"%in%args) useReleaseData=TRUE else useReleaseData=FALSE

# any special colouring? Otherwise default is genotype call.
# options for -colour are any categorical variables in the sample table.
colourOptions = c("Submitted.Gender","Inferred.Gender","Country.of.birth..UK.elsewhere.","Ethnic.background","Submitted.Plate.Name","Submitted.Well","Processed.Plate.Name","Processed.Well","Plate.Scanned.Date","Scanner")

if(useReleaseData) colourOptions =  c("Submitted.Gender","Inferred.Gender","Country.of.birth..UK.elsewhere.","Ethnic.background","Plate.Name","Submitted.Well","Plate.Name","Well","Plate.Scanned.Date","Scanner")


is.a.colourField=FALSE
if("-colour" %in% args){
    colourField = args[which(args=="-colour")+1]
    if(!colourField%in%colourOptions) print(paste0("Your choice of -colours is not an option in the sample table. Pick one of these: ",paste0(colourOptions,collapse=", ") ) ) else is.a.colourField=T
}

is.samples=FALSE
if("-samples" %in% args){
    # comma-separated list of samples IDs, or a file with them in it.
    selectSamples = args[which(args=="-samples")+1]
    samples = read.table(selectSamples,stringsAsFactors=F,header=F)[,1]
    is.samples=TRUE
    print(head(samples))
    print( paste0(length(samples)," samples will be printed within their batch." ) ) 
}

# do we separate sexes?
sexSep=FALSE
if("-sex" %in% args) sexSep=TRUE

# do we plot original A and B intensities, or the transformed ones? Default = transformed
orig.intensity=FALSE
if("-orig"%in%args) orig.intensity=TRUE


# do we plot all batches on one page???
onepage=FALSE
if("-onepage"%in%args) onepage=TRUE

# do we only plot the legend for these snps?
plotlLegendOnly=FALSE
if("-legendonly"%in%args) plotlLegendOnly=TRUE

# generate a random number for this set of plots (or use existing one)
if("-rNumber"%in%args) rNumber = as.numeric(args[which(args=="-rNumber")+1]) else rNumber = round(runif(n=1,min=0,max=100000))

# any extra titles?
if("-title"%in%args) rNumber = paste0(rNumber,"-",args[which(args=="-title")+1])

# do we plot pdf or png?
plotPDF=FALSE
if("-pdf"%in%args) plotPDF=TRUE

# do we show the counts of samples in each genotype?
showCounts = TRUE
if("-nocounts"%in%args) showCounts=FALSE

# do we show the maf of samples in each genotype (this overrides showCounts)?
showMaf = 0
if("-maf"%in%args) showMaf=1
if("-mac"%in%args) showMaf=2

# do we just show mac?
showMac = FALSE
if("-mac"%in%args) showMac=TRUE

# do we fix the number of plotted columns? (only works with -onepage)
maxColumns=5
if("-maxcol"%in%args) maxColumns=as.numeric(args[which(args=="-maxcol")+1])

# do we plot the SSPs?
plotSSP = TRUE
if("-nossp"%in%args) plotSSP = FALSE

# do we scale the image in any way?
imagescale=1
if("-scale"%in%args) imagescale = as.numeric(args[which(args=="-scale")+1])


################################################

system(paste0('mkdir ',outputdir))

nBatches = length(batches)
print(paste0("Plotting ",nBatches," batches. Is this what you expect?"))


#1. Get marker positions in binary files.

print('getting position of markers(s) in binary files...')

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



                                        # SNPs= c("rs3131962","rs12562034","Affx-52336937")

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

                                        #2. Get SSP information

sspUKBB = read.delim(paste0(inputdir,"/V2_QCed.SNPspecificPriors/UKBiobank_SSPs_source_batches.txt.gz"),sep="\t",header=T,stringsAsFactors=F,colClasses="character")
sspUKBL = read.delim(paste0(inputdir,"/V2_QCed.SNPspecificPriors/UKBiLEVEAX_SSPs_source_batches.txt.gz"),sep="\t",header=T,stringsAsFactors=F,colClasses="character")


                                        #3.
############################################
                                        # if we are using the released data, then load the sample files
if( ( useReleaseData )&(!plotlLegendOnly )){

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
    
}
############################################


#4. Plot these probesets
sexColours = c("blue","red","purple")
names(sexColours) = c("M","F","")

plotThisSnpBasic <- function(i,sexSep=FALSE,snp.by.page=TRUE,useReleaseData=FALSE,onepage=FALSE,maxColumns=5,...){

    
    snpName = allSnpInfo[i,"AffySNPID"]
    pname = allSnpInfo[i,"ProbeSetID"]
    rsID = allSnpInfo[i,"dbSNPRSID"]
    chrom = allSnpInfo[i,"Chromosome"]
    position = as.numeric(allSnpInfo[i,"Position"])
    isPAR = sum(as.numeric(allSnpInfo[i,c("IsInPAR1","IsInPAR2")]))
    bim = allSnpInfo[i,"bim"]

    print(snpName)
    
    if(chrom%in%c("MT","X","Y","XY")) is.autosomal = FALSE else is.autosomal = TRUE
    
    
    # get SSP batch UKBiobank
    sspBatchUKBB = sspUKBB$batch_id[sspUKBB$probeset_id==pname]
    
    # get SSP batch UKBiLEVE 
    sspBatchUKBL = sspUKBL$batch_id[sspUKBL$probeset_id==pname]
    
    if(plotSSP) batchesThisSnp = unique(c(batches,sspBatchUKBL,sspBatchUKBB)) else batchesThisSnp = unique(batches)
    
    nBatches = length(batchesThisSnp)

    # if it's a X-chromosome, and not in the PAR, then automatically plot males and females separately.

    if((chrom=="X") & (isPAR==0)) sexSep=T
    
    if(nBatches >= 7) arrangement = c(3,3)
    if(nBatches%in%c(5,6)) arrangement = c(2,3)
    if(nBatches%in%c(3,4)) arrangement = c(2,2)
    if(nBatches==2) arrangement = c(1,2)
    if(nBatches==1) arrangement = c(1,1)
    if(sexSep) arrangement = c(min(3,nBatches),3)
    if((is.a.colourField)&(!sexSep)) arrangement = c(min(3,nBatches),2)
    if(is.samples) arrangement = c(min(3,nBatches),1+sum(c(is.samples,is.a.colourField)))
    if(is.samples&orig.intensity) arrangement = c(min(3,nBatches),3) # ad-hoc for plotting original intensities
    scale=1
    if(onepage) {
                                        # all on onepage overrides all other arrangements.
        nImages = sum(nBatches*c(1,sexSep,is.a.colourField,is.samples,orig.intensity))
        print(paste0("printing ",nImages," plots."))
        # split over multiple pages when more than 5 columns
        nColumns = min( ceiling(sqrt(nImages)), maxColumns)
        arrangement = c(min( ceiling(nImages/nColumns),maxColumns) , nColumns)        
        print(paste0("printing ",prod(arrangement)," images on each page."))

        scale=max(arrangement[1],4)/4 # only scale if more than 4 in column
    }

    ratio = arrangement[1]/arrangement[2]
    
    qc = "qcdata-"

    if(useReleaseData){
                                        # if using the released data, read in data for each snp across all batches and subset by batch later.
                                        # chromosome names are numeric in the files, but alphanumeric in the file names.
        pid.bim = as.numeric(allSnpInfo[i,"line"])

        if(is.na(pid.bim)) {
            
            print("This snp was not found in the release data! Skipping...")
            return(NULL)
        }
        
        bimChrom = str_split(str_split(basename(bim),"chr")[[1]][2],"\\.")[[1]][1]
        
        calls = unpack.calls2(paste0(callsFiles,".chr",bimChrom,".bed"),famFile=NULL,pid.bim)
        intensities = unpack.intensities2(paste0(intensityFiles,".chr",bimChrom,".bin"),famFile=NULL,pid.bim)
        dataAll = list(calls=calls,intensities=intensities)
        
        print( paste0(callsFiles,".chr",bimChrom,".bed") )
        print( paste0(intensityFiles,".chr",bimChrom,".bin") )
        print(rsID)

        qc="reldata-"
    }

    if(plotPDF){
        pdf(file=paste0(outputdir,"/clusterplot-",qc,rNumber,"__",rsID,"__",snpName,"__",pname,"__chr",chrom,"-",position,"__",batchesPlotted,".pdf"),width=7*scale*imagescale,height=7*scale*imagescale*ratio) 
    } else {
        png(file=paste0(outputdir,"/clusterplot-",qc,rNumber,"__",rsID,"__",snpName,"__",pname,"__chr",chrom,"-",position,"__",batchesPlotted,"__%02d.png"),width=7*scale*imagescale,height=7*scale*imagescale*ratio,units="in",res=200)    
    }
    print(imagescale)
        par(mfrow=arrangement,cex=imagescale/2)
    #    par(mfrow=arrangement,cex=imagescale/1.7) #### <=== TEMP for sex!
    #   par(mfrow=arrangement,cex=imagescale/1.2) #### <=== TEMP for plate!

    marPar = par()$mar
    if(onepage) par(mfrow=arrangement,mar=c(2,2,marPar[3],marPar[4]),oma=c(3,3,0,0),cex.main=(1/scale))

    imageCount = 0
    for(Batch in sort(batchesThisSnp)){
        
        print(Batch)
        
        if( !useReleaseData ){
                                        # use original qc data
            
            if(sub.is.ukbiobank(Batch)) {
                pid = as.numeric(allSnpInfo[i,"pid.UKBB"])
            } else {
                pid = as.numeric(allSnpInfo[i,"pid.UKBL"])
            }
            
            if(is.na(pid)) {
                print(paste0("WARNING: This SNP is not available in the array for this batch: ",Batch))
                next
            }
            
            phenofile = paste0(inputdir,"/",Batch,"/V2_QCed_",Batch,"_Sample_Table_Pheno.csv")        
            datadir = paste0(inputdir,"/",Batch)
            data = unpack.probeset.data(datadir,pid=pid,is.autosomal=is.autosomal,pname=pname)

            batchInfo = get.batch.info(Batch,phenofile)                
        }

        if( useReleaseData ){
                                        # use released data sets
            
            if(sub.is.ukbiobank(Batch)) {
                pid = as.numeric(allSnpInfo[i,"pid.UKBB"])
            } else {
                pid = as.numeric(allSnpInfo[i,"pid.UKBL"])
            }
            
            if(is.na(pid)) {
                print(paste0("WARNING: This SNP is not available in the array for this batch: ",Batch))
                next
            }


            posteriors = unpack.snp.posterior2(paste0(inputdir,"/",Batch),pid,pname=pname,is.autosomal)
            
            subset = (outTable$Batch==Batch)&(!is.na(outTable$Batch))
            batchInfo = outTable[subset,] # only plot samples in this batch. outTable should be in the order of famOrder.
            intensitySubset = sort(c(which(subset)*2,which(subset)*2-1))
            data = list(calls=dataAll$calls[subset],intensities=dataAll$intensities[intensitySubset],posteriors=posteriors$posteriors,posteriors1=posteriors$posteriors1)
            
        }
                                        # is this the SSP batch?
        if(plotSSP) is.ssp = Batch%in%c(sspBatchUKBL,sspBatchUKBB) else is.ssp=FALSE
        if(!is.ssp) extra = c("","")
        if(is.ssp) extra = c("===>  ","  <===")
        
        if(length(data) == 0) print(paste0("WARNING: problem with accessing binary data for this snp/probeset ",snpName," ",pname)) else {
            
            cluster.plot(data,main.text=paste0(extra[1],Batch,extra[2]),...) # all samples together by genotype call
            imageCount=imageCount+1
            
            colourSamples=c()
            if(is.a.colourField) {
                
                colourFieldValues = batchInfo[,colourField]
                myColours = c(sample(colors(),length(unique(colourFieldValues))))
                names(myColours) = unique(colourFieldValues)

                if(colourField=="Ethnic.background") colourSamples = ethnicity2col[colourFieldValues]
                if(grepl("Gender",colourField)) colourSamples = sexColours[colourFieldValues]                
                if(!length(colourSamples)) colourSamples = myColours[colourFieldValues]
                
                if(!sexSep) {
                    cluster.plot(data,col=colourSamples,main.text=paste0(extra[1],Batch," by ",colourField,extra[2]),Alpha=1,...) # all samples together coloured by requested field.
                    imageCount=imageCount+1
                }
            }

            
            if(is.samples){
                
                colourSamples=rep("transparent",nrow(batchInfo))
                samplesToPlot = batchInfo$PIID%in%samples

                
                if(sum(samplesToPlot) > 0){
                    print("printing some samples for this batch")
                    colourSamples[samplesToPlot]="red"
                    cluster.plot(data,col=colourSamples,main.text=paste0(extra[1],Batch," ",sum(samplesToPlot)," samples",extra[2]),...)
                    imageCount=imageCount+1
                    
                } else {
                    plot.new();
                }
            }

                                        # one plot with the transformed scale.
            if(orig.intensity) {
                cluster.plot(data,main.text=paste0(extra[1],Batch,extra[2]),orig.scale=FALSE,...) # all samples together by genotype call
                imageCount=imageCount+1
            }

            
            if(sexSep) {
                print("sex")
                Gender = batchInfo$Inferred.Gender                
                cluster.plot(data,exclude.indivs=which(Gender=="F"),col=colourSamples,main.text=paste0("males"),...) # males
                imageCount=imageCount+1
                cluster.plot(data,exclude.indivs=which(Gender=="M"),col=colourSamples,main.text=paste0("females"),...) # females
                imageCount=imageCount+1
            }
        }

        if( ( onepage ) ){
            if( ( imageCount==nImages )|( imageCount==prod(arrangement)) ){
                                        # plot axis labels on outer margins of combined plots.
                print("axis labels")
                print(imageCount)
                mtext(text="Contrast:  log2(A/B)",side=1,line=1,outer=TRUE)
                mtext(text="Strength:  log2(A*B)/2",side=2,line=1,outer=TRUE)        
            }
        }
        
    }
    
    
    dev.off()
}




############################################
# 5. Plot a legend for each SNP
############################################

nsnps = nrow(allSnpInfo)

legendPlot <- function(i){

    snpName = allSnpInfo[i,"AffySNPID"]
    pname = allSnpInfo[i,"ProbeSetID"]
    rsID = allSnpInfo[i,"dbSNPRSID"]
    chrom = allSnpInfo[i,"Chromosome"]
    position = as.numeric(allSnpInfo[i,"Position"])
    isPAR = sum(as.numeric(allSnpInfo[i,c("IsInPAR1","IsInPAR2")]))

#    par(cex=3)
    
    plot(NULL,xlim=c(0,1),ylim=c(0,1),axes=FALSE,xlab=NA,ylab=NA)
    if(useReleaseData) alleles = allSnpInfo[i,c("Reference","Alternative")] else alleles = allSnpInfo[i,c("AlleleA","AlleleB")]
    
    if( chrom%in%c("MT","Y") ) {
        leg.call = list("legend"=c("A","B","no call","",paste0(":= ",alleles)),"col"=c(geno.cols[c(2,4,1)],NA,"black","black"),"pch"=c(19,19,4,NA,65,66))
    } else {
        leg.call = list("legend"=c("AA","AB","BB","no call","",paste0(":= ",alleles)),"col"=c(geno.cols[c(2,3,4,1)],NA,"black","black"),"pch"=c(19,19,19,4,NA,65,66))
    }
    
    legend("topleft",legend=leg.call$legend,col=leg.call$col,pch=leg.call$pch,bty="n",xpd=NA,title=paste0("Marker: ", rsID,"\nPosition: chr",chrom,":",position),title.adj=0)
    
    return(NULL)
}

print("printing legend for these snps/plots...")

if(plotPDF){
    pdf(file=paste0(outputdir,"/clusterplot-",rNumber,"-LEGEND.pdf"),width=7,height=7,bg="transparent")
} else {
    png(file=paste0(outputdir,"/clusterplot-",rNumber,"-LEGEND.%02d.png"),width=7,height=7,units="in",res=200,bg="transparent")
}

par(mfrow=c(5,5))
sapply(1:nsnps,legendPlot)

dev.off()


##########################################
# Actually call the plotting function

if(!plotlLegendOnly) {


    
    if(nsnps > 50) {
                                        #print("You have asked to plot more than 50 SNPs. Instead I'll plot the first 50")
                                        #nsnps=50
    }

    if(!onepage) out = sapply(1:nsnps,function(x) plotThisSnpBasic(x,orig.scale=orig.intensity,sexSep=sexSep,useReleaseData=useReleaseData,onepage=FALSE,showCounts=showCounts,showMaf=showMaf))

    if(onepage) out = sapply(1:nsnps,function(x) plotThisSnpBasic(x,orig.scale=orig.intensity,sexSep=sexSep,useReleaseData=useReleaseData,onepage=TRUE,xlab.text=NA,ylab.text=NA,showCounts=showCounts,showMaf=showMaf,maxColumns=maxColumns))
    
}


print(paste0("DONE! plots saved in ",outputdir,"/clusterplots-",rNumber,"*"))
