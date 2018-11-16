## Script to analyse PCA with respect to ethnic background or country of birth

library(zoo)
args=commandArgs(trailingOnly=T)

#args = c('/well/ukbiobank/expt/V2_QCed.SNP-QC/src/V2_QCed.snpqc-tests.R','/well/ukbiobank/expt/V2_QCed.SNP-QC/src/V2_QCed.bin2clusterplots.R','/well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/R/scripts/readPSperformance.R','/well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/R/scripts/auxFunctions.R','-in','b1__b11-b001__b095-pca-UKbio-round2')


print(args)
h = args[-c(which(args%in%c("-in")),1+which(args%in%c("-in")))]
for(helperScript in h){
    source(helperScript)
}


#################  SnpLoads should be the same file used in the submit-pca-UKBio.sh script.
BatchesFile = paste0(baseSampleQCDir,"/QC-Scripts/batchList.txt")
inputFile = args[which(args%in%c("-in"))+1]
OutDir =  paste0(baseSampleQCDir,"/QC-Scripts/PCA/pca-UKBio")
Rfile=paste0(baseSampleQCDir,"/data/PCA/",inputFile,".RData")
#################

system(paste0('mkdir ',OutDir,'/plots'))

BatchList = read.table(BatchesFile,header = FALSE, stringsAsFactors = FALSE)
nBatches = nrow(BatchList)

Info=read.multiple.batch.info(batchInfoFields[c(45,43,36,55:59)])
Info = dplyr::tbl_df(Info)

load(Rfile,verbose=TRUE)
load(paste0(OutDir,"/pca-UKBio-countryColours.RData"),verbose=TRUE)

Info = left_join(PCs[,"PIID"],Info,by=c("PIID"="PIID"))



#################
# get country of birth and ethnicity information

Info$Country.of.birth=get.place.of.birth(Info)
Info$Eth2 = ethnicity2pop(Info$Ethnic.background)
cob = Info$Country.of.birth
eth2 = Info$Eth2
cob2 = cob
tab= table(cob)
minSamples=20
tooSmall = names(tab[tab<minSamples])
cob2[cob%in%tooSmall] = NA

print(paste0(sum(cob%in%tooSmall)," individuals removed because there is fewer than ",minSamples," individuals with their birth country. ",100*sum(cob%in%tooSmall)/length(cob),"% of data."))

nPCs = sum(grepl("^PC",colnames(PCs)))


#################
# Apply linear model for each PC with ethnic background and country of birth covariates


print(paste0("Running linear regression on ",nPCs," PCs for ethnic background..."))

ethLMnorm = sapply(1:nPCs,function(pc){
    print(pc)
    y = PCs[[paste0("PC",pc)]]
    yNorm = y/sd(y)  # normalise by standard deviation and mean
#    LMfit = lm(y~0+eth2+cob)
    LMfit = lm(yNorm~0+eth2)   
    s = summary(LMfit)
    return(s)
},simplify=FALSE)

ethLMraw = sapply(1:nPCs,function(pc){
    print(pc)
    y = PCs[[paste0("PC",pc)]]
    #yNorm = y/sd(y)  # normalise by standard deviation and mean
#    LMfit = lm(y~0+eth2+cob)
    LMfit = lm(y~0+eth2)   
    s = summary(LMfit)
    return(s)
},simplify=FALSE)

# Allow the intercept to change each time (although will probably always be the mean of British people!)
ethLMnormInt = sapply(1:nPCs,function(pc){
    print(pc)
    y = PCs[[paste0("PC",pc)]]
    yNorm = y/sd(y)  # normalise by standard deviation and mean
#    LMfit = lm(y~0+eth2+cob)
    LMfit = lm(yNorm~eth2)   
    s = summary(LMfit)
    return(s)
},simplify=FALSE)


signif = 0.05/prod(dim(coeffMat))

# normalised coefficient
coeffMat = sapply(ethLMnorm,function(s) s$coefficients[,1])
pvals = sapply(ethLMnorm,function(s) s$coefficients[,4])
coeffMat2 = coeffMat
coeffMat2[pvals>=signif] = 0

# raw coefficients
coeffMatraw = sapply(ethLMraw,function(s) s$coefficients[,1])
pvalsRaw = sapply(ethLMraw,function(s) s$coefficients[,4])
coeffMatraw2 = coeffMatraw
coeffMatraw2[pvalsRaw>=signif] = 0

# normalised no intercept coefficients
coeffMatInt = sapply(ethLMnormInt,function(s) s$coefficients[,1])
pvalsInt = sapply(ethLMnormInt,function(s) s$coefficients[,4])
coeffMatInt2 = coeffMatInt
coeffMatInt2[pvalsInt>=signif] = 0


rownames(coeffMat) = rownames(coeffMat2) = rownames(pvals) = rownames(coeffMatraw) = rownames(coeffMatraw2) = rownames(pvalsRaw) = rownames(coeffMatInt) = rownames(coeffMatInt2) = rownames(pvalsInt) = gsub("eth2","",rownames(coeffMat))
colnames(coeffMat) = colnames(coeffMat2) = colnames(pvals) = colnames(coeffMatInt) = colnames(coeffMatInt2) = colnames(pvalsInt) = colnames(coeffMatraw) = colnames(coeffMatraw2) = colnames(pvalsRaw) = paste0("PC",1:nPCs)

#test = lm(as.matrix(pcs)~eth2)


#################
# Plot the matrix of coefficients

Legend = cbind.data.frame(names(ethnicity2col),ethnicity2col,ethnicity2char[names(ethnicity2col)])  # make legend
colnames(Legend) = c("Factor","colour","shape")
rownames(Legend) = names(ethnicity2col)
    
treeOrder = hclust(dist(coeffMatraw)) # order by raw coefficients

toPlot = t(coeffMat2[treeOrder$order,])
capValue = mean(abs(toPlot[toPlot!=0])) + 2*sd(abs(toPlot[toPlot!=0]))

plotMixtureHeat(toPlot,filename=paste0(OutDir,"/plots/",inputFile,"-ethnicBackgroundEffectNormalised.png"),cap=T,capValue=capValue,
title="PCy ~ 0 + ethnicty-x",xlab="PC",ylab=NA,colLabels=colnames(toPlot),rowLabels=1:nPCs,scaleTitle="Effect size (per sd.)",Legend=Legend,plotPoints=T,
height=1200,width=2000,res=150)


toPlot = t(coeffMatraw2[treeOrder$order,])
capValue = mean(abs(toPlot[toPlot!=0])) + 2*sd(abs(toPlot[toPlot!=0]))

plotMixtureHeat(toPlot,filename=paste0(OutDir,"/plots/",inputFile,"-ethnicBackgroundEffectRaw.png"),cap=T,capValue=capValue,
title="PCy ~ 0 + ethnicty-x",xlab="PC",ylab=NA,colLabels=colnames(toPlot),rowLabels=1:nPCs,scaleTitle="Effect size",Legend=Legend,plotPoints=T,
height=1200,width=2000,res=150)


toPlot = t(coeffMatInt2[treeOrder$order,])
capValue = mean(abs(toPlot[toPlot!=0])) + 2*sd(abs(toPlot[toPlot!=0]))

plotMixtureHeat(toPlot,filename=paste0(OutDir,"/plots/",inputFile,"-ethnicBackgroundEffectNormalisedIntercept.png"),cap=T,capValue=capValue,
title="PCy ~ 0 + ethnicty-x",xlab="PC",ylab=NA,colLabels=colnames(toPlot),rowLabels=1:nPCs,scaleTitle="Effect size",Legend=Legend,plotPoints=T,
height=1200,width=2000,res=150)






#################
# Apply linear model for each PC with country of birth as covariates

print(paste0("Running linear regression on ",nPCs," PCs for country of birth..."))

cobLMnorm = sapply(1:nPCs,function(pc){
    print(pc)
    y = PCs[[paste0("PC",pc)]]
    yNorm = y/sd(y)  # normalise by standard deviation
                                        #    LMfit = lm(y~0+eth2+cob)
    LMfit = lm(yNorm~0+cob2)   
    s = summary(LMfit)
    return(s)
},simplify=FALSE)


cobLMraw = sapply(1:nPCs,function(pc){
    print(pc)
    y = PCs[[paste0("PC",pc)]]
                                        #yNorm = y/sd(y)  # normalise by standard deviation and mean
                                        #    LMfit = lm(y~0+eth2+cob)
    LMfit = lm(y~0+cob2)   
    s = summary(LMfit)
    return(s)
},simplify=FALSE)

uCob2 = unique(cob2);uCob2 = uCob2[!is.na(uCob2)] 
cobLMnormMean1 = sapply(1:nPCs,function(pc){
    print(pc)
    y = PCs[[paste0("PC",pc)]]
    yNorm = y/sd(y)  # normalise by standard deviation
    o = sapply(uCob2,function(u) mean(yNorm[cob2==u],na.rm=TRUE))
    return(o)
},simplify=TRUE)


coeffMat = sapply(cobLMnorm,function(s) s$coefficients[,1])
pvals = sapply(cobLMnorm,function(s) s$coefficients[,4])

signif = 0.05/prod(dim(coeffMat))
coeffMat2 = coeffMat
coeffMat2[pvals>=signif] = 0

                                        # raw coefficients
coeffMatraw = sapply(cobLMraw,function(s) s$coefficients[,1])
pvalsRaw = sapply(cobLMraw,function(s) s$coefficients[,4])
coeffMatraw2 = coeffMatraw
coeffMatraw2[pvalsRaw>=signif] = 0

                                        # just the means
rownames(cobLMnormMean1) = uCob2
cobLMnormMean = cobLMnormMean1[match(gsub("cob2","",rownames(coeffMat)),rownames(cobLMnormMean1)),]
cobLMnormMean2 = cobLMnormMean
cobLMnormMean2[pvals>=signif] = 0


rownames(coeffMat) = rownames(coeffMat2) = rownames(pvals) = rownames(coeffMatraw) = rownames(coeffMatraw2) = rownames(pvalsRaw) = gsub("cob2","",rownames(coeffMat))
colnames(coeffMat) = colnames(coeffMat2) = colnames(pvals) = colnames(coeffMatraw) = colnames(pvalsRaw) = paste0("PC",1:nPCs)


#################
                                        # Plot the matrix of coefficients
cobCols = getColoursDistant3(length(table(cob2)))
names(cobCols) = names(sort(table(cob2)))
cobChars = rep(c(1:5),100)[1:length(cobCols)]
names(cobChars) = names(cobCols)

                                        # Make colours based on most common ethnicity value
maxEthnic = t(sapply(names(sort(table(cob2))), function(country){
    if(country=="Other/Unknown") return(list("Other/Unknown",1))
    tab = table(eth2[cob2==country])
    tab = tab[!is.na(tab)]
                                        #tab = tab[names(tab)!="Other/Unknown"]
    if(length(tab)==0) print(country)
    if(length(tab)>1) dist = (max(tab) - max(tab[tab!=max(tab)]))/sum(tab) else dist = 1
                                        #if(dist < 0.5) print(country)
    maxEth = names(tab)[tab==max(tab)][1]
                                        # }
    return(list(maxEth,dist))
},simplify=TRUE))

cobCols = ethnicity2col[unlist(maxEthnic[,1])]
cobChars = ethnicity2char[unlist(maxEthnic[,1])]
names(cobCols) = names(cobChars) = rownames(maxEthnic)


Legend = cbind.data.frame(unlist(maxEthnic[,1]),cobCols,cobChars[names(cobCols)])  # make legend
colnames(Legend) = c("Factor","colour","shape")
rownames(Legend) = names(cobCols)


treeOrder = hclust(dist(coeffMat)) # order by normalised coefficients
treeOrder = hclust(dist(coeffMat2)) # order by normalised coefficients after exluding non-significant elements


###########
                                        # Save these results!
###########

save(cobLMnormMean2,cobLMnormMean,coeffMat2,Legend,maxEthnic,signif,pvals,treeOrder,file=paste0(OutDir,"/",inputFile,"-COBbyPCs-lm.RData"))


###########
# Plot them
###########

toPlot = t(coeffMat2[treeOrder$order,])
capValue = mean(abs(toPlot[toPlot!=0])) + 2*sd(abs(toPlot[toPlot!=0]))

plotMixtureHeat(toPlot,filename=paste0(OutDir,"/plots/",inputFile,"-COBEffectNormalised.png"),cap=T,capValue=capValue,
title="PCy ~ 0 + countryOfBirth-x",xlab="PC",ylab=NA,colLabels=colnames(toPlot),rowLabels=1:nPCs,scaleTitle="Effect size (per sd.)",Legend=Legend,plotPoints=T,
height=2000,width=1500,res=150,Cex=0.9,y.cex=0.7,ylabMar=10,scaleWidth=0.1)



toPlot = t(coeffMatraw[treeOrder$order,])
capValue = mean(abs(toPlot[toPlot!=0])) + 2*sd(abs(toPlot[toPlot!=0]))

plotMixtureHeat(toPlot,filename=paste0(OutDir,"/plots/",inputFile,"-COBEffectRaw.png"),cap=T,capValue=capValue,
title="PCy ~ 0 + countryOfBirth-x",xlab="PC",ylab=NA,colLabels=colnames(toPlot),rowLabels=1:nPCs,scaleTitle="Effect size",Legend=Legend,plotPoints=T,
height=2000,width=1500,res=150,Cex=0.9,y.cex=0.7,ylabMar=10,scaleWidth=0.1)


toPlot = t(cobLMnormMean2[treeOrder$order,])
capValue = mean(abs(toPlot[toPlot!=0])) + 2*sd(abs(toPlot[toPlot!=0]))

plotMixtureHeat(toPlot,filename=paste0(OutDir,"/plots/",inputFile,"-COBNormalisedMeans.png"),cap=T,capValue=capValue,
title="PCy ~ 0 + countryOfBirth-x",xlab="PC",ylab=NA,colLabels=colnames(toPlot),rowLabels=1:nPCs,scaleTitle="sd-scaled mean",Legend=Legend,plotPoints=T,
height=2000,width=1500,res=150,Cex=0.9,y.cex=0.7,ylabMar=10,scaleWidth=0.1)
