########
# usage:
# Rscript /well/ukbiobank/expt/V2_QCed.SNP-QC/src/V2_QCed.makeClusterPlots.R <input directory with batch-based binary files> <comma-separated list of Affy SNPIDs, or file name with list> <output directory> [-b <batchName, or range, or comma-separated list>] [-colour <variable in sample table. See colourOptions>] [-orig]
#
# DEFAULTS: prints snp at all batches, 9 batches to a page
########

library(stringr)

source('/well/ukbiobank/expt/V2_QCed.SNP-QC/src/V2_QCed.bin2clusterplots.R')
source('/well/ukbiobank/expt/V2_QCed.SNP-QC/src/V2_QCed.batch2sampleinfo.R')

args = commandArgs(trailingOnly=T)
#args = c("/well/ukbiobank/expt/V2_QCed.SNP-QC/data","Affx-13890966,Affx-13890977,Affx-13890978,Affx-52187954,Affx-80256891","/well/ukbiobank/expt/V2_QCed.SNP-QC/data/plots/test_clusterPlots_plate","-b","Batch_b080-Batch_b085","-colour","Processed.Plate.Name","-samples","/well/ukbiobank/qcoutput.V2_QCed.sample-QC-testing/QC-Scripts/PCA/pca-UKBio/b1__b11-b001__b095-pca-UKbio-proj40-dim50-pc17-snpsWithHighloads0.05.txt")

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


# any special colouring? Otherwise default is genotype call.
# options for -colour are any categorical variables in the sample table.
colourOptions = c("Submitted.Gender","Inferred.Gender","Country.of.birth..UK.elsewhere.","Ethnic.background","Submitted.Plate.Name","Submitted.Well","Processed.Plate.Name","Processed.Well","Plate.Scanned.Date","Scanner")

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

########

system(paste0('mkdir ',outputdir))

nBatches = length(batches)
print(paste0("Plotting ",nBatches," batches. Is this what you expect?"))


#1. Get marker positions in binary files.

print('getting position of markers(s) in binary files...')

getSnpInfo <- function(SNPs,ps2snpTableAutosome,ps2snpTableSexchrom){

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

                                        # check that all SNPs are represented here
    if(sum(!SNPs%in%allSnpInfo[,"AffySNPID"]) > 0 ) {
        
        ids = SNPs[!SNPs%in%allSnpInfo[,"AffySNPID"]]
        print(paste0("WARNING: Unable to find markers ",paste(ids,collapse=", ")," in the axiom chip tables ",ps2snpTableAutosome," or ",ps2snpTableSexchrom,". Check that you have specified the right chip file. Will obviously not plot these markers (in some batches)."))
        
    }

    return(allSnpInfo)
}


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



#3. Plot these probesets
sexColours = c("blue","red","purple")
names(sexColours) = c("M","F","")

plotThisSnpBasic <- function(i,sexSep=FALSE,snp.by.page=TRUE,...){
    
    snpName = allSnpInfo[i,"AffySNPID"]
    pname = allSnpInfo[i,"ProbeSetID"]
    rsID = allSnpInfo[i,"dbSNPRSID"]
    chrom = allSnpInfo[i,"Chromosome"]
    position = as.numeric(allSnpInfo[i,"Position"])
    isPAR = sum(as.numeric(allSnpInfo[i,c("IsInPAR1","IsInPAR2")]))
    
    print(snpName)
    
    if(chrom%in%c("MT","X","Y")) is.autosomal = FALSE else is.autosomal = TRUE
    
    
    # get SSP batch UKBiobank
    sspBatchUKBB = sspUKBB$batch_id[sspUKBB$probeset_id==pname]
    
    # get SSP batch UKBiLEVE 
    sspBatchUKBL = sspUKBL$batch_id[sspUKBL$probeset_id==pname]
    
    batchesThisSnp = unique(c(batches,sspBatchUKBL,sspBatchUKBB))
    
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

    ratio = arrangement[1]/arrangement[2]
    
    png(file=paste0(outputdir,"/clusterplot__",snpName,"__",pname,"__chr",chrom,"-",position,"__",batchesPlotted,"__%02d.png"),width=7,height=7*ratio,units="in",res=150)    

    par(mfrow=arrangement)
    
    for(Batch in batchesThisSnp){        
      
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
       
        # is this the SSP batch?
        is.ssp = Batch%in%c(sspBatchUKBL,sspBatchUKBB)
        if(!is.ssp) extra = c("","")
        if(is.ssp) extra = c("===>  ","  <===")
        
        if(length(data) == 0) print(paste0("WARNING: problem with accessing binary data for this snp/probeset ",snpName," ",pname)) else {
            
            cluster.plot(data,main.text=paste0(extra[1],Batch,extra[2]),...) # all samples together by genotype call
            colourSamples=c()
            if(is.a.colourField) {
                
                colourFieldValues = batchInfo[,colourField]
                myColours = c(sample(colors(),length(unique(colourFieldValues))))
                names(myColours) = unique(colourFieldValues)

                if(colourField=="Ethnic.background") colourSamples = ethnicity2col[colourFieldValues]
                if(grepl("Gender",colourField)) colourSamples = sexColours[colourFieldValues]                
                if(!length(colourSamples)) colourSamples = myColours[colourFieldValues]
                
                if(!sexSep) cluster.plot(data,col=colourSamples,main.text=paste0(extra[1],Batch," by ",colourField,extra[2]),...) # all samples together coloured by requested field.                
            }

            
            if(is.samples){
                
                colourSamples=rep("transparent",nrow(batchInfo))
                samplesToPlot = batchInfo$PIID%in%samples

                
                if(sum(samplesToPlot) > 0){
                    print("printing some samples for this batch")
                    colourSamples[samplesToPlot]="red"
                    cluster.plot(data,col=colourSamples,main.text=paste0(extra[1],Batch," ",sum(samplesToPlot)," samples",extra[2]),...)
                } else {
                    plot.new();
                }
            }

                                                    # one plot with the transformed scale.
            if(orig.intensity) cluster.plot(data,main.text=paste0(extra[1],Batch,extra[2]),orig.scale=FALSE) # all samples together by genotype call

            
            if(sexSep) {
                print("sex")
                Gender = batchInfo$Inferred.Gender                
                cluster.plot(data,exclude.indivs=which(Gender=="F"),col=colourSamples,main.text=paste0("males"),...) # males
                cluster.plot(data,exclude.indivs=which(Gender=="M"),col=colourSamples,main.text=paste0("females"),...) # females
            }
        }
    }        
    dev.off()
}


nsnps = nrow(allSnpInfo)

if(nsnps > 50) {
    #print("You have asked to plot more than 50 SNPs. Instead I'll plot the first 50")
    #nsnps=50
}

out = sapply(1:nsnps,function(x) plotThisSnpBasic(x,orig.scale=orig.intensity,sexSep=sexSep))


print(paste0("DONE! plots saved in ",outputdir))
