## Script to compute PC-adjusted heterozygosity (and missingness)

library(dplyr)

args = commandArgs(trailingOnly=T)
#args = c("/well/ukbiobank/expt/V2_QCed.SNP-QC/src/V2_QCed.snpqc-tests.R","/well/ukbiobank/expt/V2_QCed.SNP-QC/src/V2_QCed.bin2clusterplots.R","/well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/R/scripts/readPSperformance.R","/well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/R/scripts/auxFunctions.R",
#"-in","/well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/Combined/b1__b11-b001__b095-autosome-sampleqc","-pcs","/well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/PCA/b1__b11-b001__b095-pca-UKbio-init.RData","-npcs","6PCs")


print(args)
h = args[-c(which(args%in%c("-in","-pcs","-npcs")),1+which(args%in%c("-in","-pcs","-npcs")))]
for(helperScript in h){
    source(helperScript)
}

input = args[which(args=="-in") + 1]
outname = basename(input)
nPCs = args[which(args=="-npcs")+1]

# get PCs info
pcsFile = args[which(args=="-pcs") + 1]
load(pcsFile,verbose=T)

# get het info
hetFile = paste0(input,".het.gz")
hetero = read.table(gzfile(hetFile),stringsAsFactors=FALSE,header=T)

# get missing info
missFile = paste0(input,".imiss.gz")
miss = read.table(gzfile(missFile),stringsAsFactors=FALSE,header=T)

if( sum(miss$IID!=hetero$IID) > 0 ) print("WARNING: samples in het and imiss files are not in the same order...")

# compute heterozygosity
hetero$het = (hetero$N.NM. - hetero$O.HOM.)/hetero$N.NM.
hetero$miss = miss$F_MISS
hetero$logit.miss = log( hetero$miss / (1 - hetero$miss) )
hetero = dplyr::tbl_df(hetero)

# join with PCs
Table = left_join(hetero,PCs,by=c("IID"="PIID"))

if(dim(PCs)[1] != dim(hetero)[1]) print( paste0("WARNING: Het and PCs files have different numbers of samples. Do you expect this? Only using samples in het and missing files: ",input) )

if(dim(PCs)[1] != dim(Table)[1]) print( paste0("WARNING: Het and PCs files have different numbers of samples. Do you expect this? Only using samples in het and missing files: ",input) )

    
#########################################################################
# Functions for computing het correction

compute.het.correction.4PCs = function(hetmisspcs) {

    print("Correcting using the first 4 PCs")

  ## Suppose that ancestry can be inferred (or more generally, is correlated with)
  ## the first 4 principal components of the genotype matrix
  if (!min(c("het","PC1","PC2","PC3","PC4") %in% names(hetmisspcs))) {
    return(NULL)
  }
  frmla = formula("het ~ (PC1 + PC2 + PC3 + PC4)^2 + I(PC1*PC1) + I(PC2*PC2) + I(PC3*PC3) + I(PC4*PC4)")
  lmfit = lm(formula = frmla, data = hetmisspcs)
  betahat = coefficients(lmfit)
  Intercept = betahat[1]
  return (list(frmla=frmla,lmfit=lmfit,betahat=betahat,Intercept=Intercept))
}
apply.het.correction.4PCs = function(hetmisspcs,het.correct) {
  ## Suppose that ancestry can be inferred (or more generally, is correlated with)
  ## the first 4 principal components of the genotype matrix
  if (!min(c("het","PC1","PC2","PC3","PC4") %in% names(hetmisspcs))) {
    return(NULL)
  }
  het = hetmisspcs$het
  lmfit = het.correct$lmfit
  betahat = het.correct$betahat
  Intercept = het.correct$Intercept
  ## het.corrected = het - (model.matrix(lmfit) %*% betahat - Intercept)
  het.corrected = het - (predict(lmfit, newdata = hetmisspcs) - Intercept)
  return (het.corrected)
}


# Need to add in a 5th PC to capture south American-based PC (PC5)

compute.het.correction.5PCs = function(hetmisspcs) {
    
    print("Correcting using the first 5 PCs")

  ## Suppose that ancestry can be inferred (or more generally, is correlated with)
  ## the first 4 principal components of the genotype matrix
  if (!min(c("het","PC1","PC2","PC3","PC4","PC5") %in% names(hetmisspcs))) {
    return(NULL)
  }
  frmla = formula("het ~ (PC1 + PC2 + PC3 + PC4 + PC5)^2 + I(PC1*PC1) + I(PC2*PC2) + I(PC3*PC3) + I(PC4*PC4) + I(PC5*PC5)")
  lmfit = lm(formula = frmla, data = hetmisspcs)
  betahat = coefficients(lmfit)
  Intercept = betahat[1]
  return (list(frmla=frmla,lmfit=lmfit,betahat=betahat,Intercept=Intercept))
}

apply.het.correction.5PCs = function(hetmisspcs,het.correct) {
  ## Suppose that ancestry can be inferred (or more generally, is correlated with)
  ## the first 4 principal components of the genotype matrix
  if (!min(c("het","PC1","PC2","PC3","PC4","PC5") %in% names(hetmisspcs))) {
    return(NULL)
  }
  het = hetmisspcs$het
  lmfit = het.correct$lmfit
  betahat = het.correct$betahat
  Intercept = het.correct$Intercept
  ## het.corrected = het - (model.matrix(lmfit) %*% betahat - Intercept)
  het.corrected = het - (predict(lmfit, newdata = hetmisspcs) - Intercept)
  return (het.corrected)
}


## 6 PCs!!

compute.het.correction.6PCs = function(hetmisspcs) {

    print("Correcting using the first 6 PCs")
    
  ## Suppose that ancestry can be inferred (or more generally, is correlated with)
  ## the first 6 principal components of the genotype matrix
  if (!min(c("het","PC1","PC2","PC3","PC4","PC5","PC6") %in% names(hetmisspcs))) {
      print("HELP! PCs file not formatted properly...")
    return(NULL)
  }
  frmla = formula("het ~ (PC1 + PC2 + PC3 + PC4 + PC5 + PC6)^2 + I(PC1*PC1) + I(PC2*PC2) + I(PC3*PC3) + I(PC4*PC4) + I(PC5*PC5) + I(PC6*PC6)")
  lmfit = lm(formula = frmla, data = hetmisspcs)
  betahat = coefficients(lmfit)
  Intercept = betahat[1]
  return (list(frmla=frmla,lmfit=lmfit,betahat=betahat,Intercept=Intercept))
}
apply.het.correction.6PCs = function(hetmisspcs,het.correct) {
  ## Suppose that ancestry can be inferred (or more generally, is correlated with)
  ## the first 4 principal components of the genotype matrix
  if (!min(c("het","PC1","PC2","PC3","PC4","PC5","PC6") %in% names(hetmisspcs))) {
    return(NULL)
  }
  het = hetmisspcs$het
  lmfit = het.correct$lmfit
  betahat = het.correct$betahat
  Intercept = het.correct$Intercept
  ## het.corrected = het - (model.matrix(lmfit) %*% betahat - Intercept)
  het.corrected = het - (predict(lmfit, newdata = hetmisspcs) - Intercept)
  return (het.corrected)
}

#########################################################################


## Apply the functions
## Should you fancy saving the beta coefficients, this is the fitted linear model

# don't adjust anything that has 0 heterozygosity (as in X chromosomes)
if( sum(Table$het==0) > 0 ){
    print("Some heterozygosity is zero. Excluding these from the table.")
    Table2 = filter(Table,het!=0)
    
    het.correct = do.call( paste0("compute.het.correction.",nPCs), list(Table2) )
    het.corrected = do.call( paste0("apply.het.correction.",nPCs), list(Table2,het.correct) )

    Table2$het.corrected = het.corrected
    Table = left_join(Table,select(Table2,het.corrected,IID),by=c("IID"="IID"))
    
} else {

    het.correct = do.call( paste0("compute.het.correction.",nPCs), list(Table) )
    het.corrected = do.call( paste0("apply.het.correction.",nPCs), list(Table,het.correct) )

    Table$het.corrected = het.corrected
    
}
save(het.correct, Table, file = paste0(input,"-hetcorrected-",nPCs,"-imiss.RData"))

gz1 <- gzfile(paste0(input,"-hetcorrected-",nPCs,"-imiss.txt.gz"), "w")
forOut = dplyr::select(Table,IID,het,het.corrected,miss,logit.miss)
write.table(forOut,file=gz1,quote=F,row.names=F,col.names=T)

print("DONE! Output files saved in:")
print( paste0(input,"-hetcorrected-",nPCs,"-imiss.RData") )
print( paste0(input,"-hetcorrected-",nPCs,"-imiss.txt.gz") )
