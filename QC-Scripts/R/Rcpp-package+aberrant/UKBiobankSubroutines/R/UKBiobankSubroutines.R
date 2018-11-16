



SNPHWE <- function(x) {
    pvalue = .Call("SNPHWE",as.integer(x))
    if (pvalue<0||pvalue>1) { pvalue = NA }
    return(pvalue)
}
Fisher2x2 <- function(x) {
    pvalue = .Call("Fisher2x2",as.integer(x))
    if (pvalue<0||pvalue>1) { pvalue = NA }
    return(pvalue)
}
Fisher2x2.Onesided <- function(x,x11.is.greater.alt) {
    pvalue = .Call("Fisher2x2_onesided",as.integer(x),as.integer(x11.is.greater.alt))
    if (pvalue<0||pvalue>1) { pvalue = NA }
    return(pvalue)
}
Fisher2x3 <- function(x) {
    pvalue = .Call("Fisher2x3",as.integer(x))
    if (pvalue<0||pvalue>1) { pvalue = NA }
    return(pvalue)
}

HardyWeinbergPvals <- function(data) {
    pvalues = rep(NA,3)
    if (sum(is.na(data))) { return(pvalues) }
    chrom = data[1]
    fcounts = data[2:5]
    mcounts = data[6:9]
    ## Plink uses the following codes to specify non-autosome chromosome types:
    ## X  -> 23 (X chromosome)
    ## Y  -> 24 (Y chromosome)
    ## XY -> 25 (Pseudo-autosomal region)
    ## MT -> 26 (Mitochondrial)
    if (chrom==23) {
        pvalues = .Call("TestHardyWeinbergX",as.integer(fcounts),as.integer(mcounts))
    } else if (chrom==25 || (chrom>=1 && chrom<=22)) { ## PAR region or autosome
        pvalues = .Call("TestHardyWeinberg",as.integer(fcounts),as.integer(mcounts))
    }
    pvalues[(pvalues<0)|(pvalues>1)] = NA
    return(pvalues)
}

CompareGenotypes <- function(data) {
    pvalues = rep(NA,3)
    if (sum(is.na(data))) { return(pvalues) }
    chrom = data[1]
    fcounts1 = data[ 2:5 ]
    mcounts1 = data[ 6:9 ]
    fcounts2 = data[10:13]
    mcounts2 = data[14:17]
    ## Plink uses the following codes to specify non-autosome chromosome types:
    ## X  -> 23 (X chromosome)
    ## Y  -> 24 (Y chromosome)             ## None of the tests is applicable to Y SNPs
    ## XY -> 25 (Pseudo-autosomal region)
    ## MT -> 26 (Mitochondrial)
    if (chrom==23) {
      pvalues = .Call("TestGenotypeFreqsX",as.integer(fcounts1),as.integer(mcounts1),
		      as.integer(fcounts2),as.integer(mcounts2))
    } else if (chrom==24) {
      pvalues = .Call("TestGenotypeFreqsY",as.integer(fcounts1),as.integer(mcounts1),
		      as.integer(fcounts2),as.integer(mcounts2))
    } else if (chrom==26) {
      pvalues = .Call("TestGenotypeFreqsMT",as.integer(fcounts1),as.integer(mcounts1),
		      as.integer(fcounts2),as.integer(mcounts2))
    } else { ## PAR region or autosome
      pvalues = .Call("TestGenotypeFreqs",as.integer(fcounts1),as.integer(mcounts1),
		      as.integer(fcounts2),as.integer(mcounts2))
    }
    pvalues[(pvalues<0)|(pvalues>1)] = NA
    return(pvalues)
}

genotype.counts <- function(calls) {
    if (is.data.frame(calls))                 { calls = as.matrix(calls) }
    if (is.matrix(calls)&&min(dim(calls))==1) { calls = as.vector(calls) }
    if (is.matrix(calls)) {
        counts = t(apply(calls,1,genotype.counts.vector))
        counts = data.frame(nAA=counts[,1],nAB=counts[,2],nBB=counts[,3],n..=counts[,4])
    } else {
        counts = genotype.counts.vector(calls)
        counts = data.frame(nAA=counts[1],nAB=counts[2],nBB=counts[3],n..=counts[4])
    }
    return(counts)
}
genotype.counts.vector <- function(x) {
    counts = .Call("GenoCounts",as.integer(x))
    return(counts)
}
calls.FisherExact <- function(data,gender,category) {
    category = as.character(category)
    chrom = data[1]
    calls = data[-1]
    fAA = 1*((calls==0)&(gender=="F"))
    mAA = 1*((calls==0)&(gender=="M"))
    fAB = 1*((calls==1)&(gender=="F"))
    mAB = 1*((calls==1)&(gender=="M"))
    fBB = 1*((calls==2)&(gender=="F"))
    mBB = 1*((calls==2)&(gender=="M"))
    if (chrom==23) {         ## X chromosome
        A = tapply(2*fAA+fAB+mAA,category,sum)
        B = tapply(2*fBB+fAB+mBB,category,sum)
        totA = sum(A)
        totB = sum(B)
        pval = apply(cbind(A,B,totA-A,totB-B),1,Fisher2x2)
    } else if (chrom==24) {  ## Y chromosome
        A = tapply(mAA,category,sum)
        B = tapply(mBB,category,sum)
        totA = sum(A)
        totB = sum(B)
        pval = apply(cbind(A,B,totA-A,totB-B),1,Fisher2x2)
    } else if (chrom==26) {  ## mitochondria
        A = tapply(fAA+mAA,category,sum)
        B = tapply(fBB+mBB,category,sum)
        totA = sum(A)
        totB = sum(B)
        pval = apply(cbind(A,B,totA-A,totB-B),1,Fisher2x2)
    } else { ## PAR region or autosome
        AA = tapply(fAA+mAA,category,sum)
        AB = tapply(fAB+mAB,category,sum)
        BB = tapply(fBB+mBB,category,sum)
        AA = tapply(calls==0,category,sum)
	AB = tapply(calls==1,category,sum)
	BB = tapply(calls==2,category,sum)
	totAA = sum(AA)
	totAB = sum(AB)
	totBB = sum(BB)
        pval = apply(cbind(AA,AB,BB,totAA-AA,totAB-AB,totBB-BB),1,Fisher2x3)
    }
    return(pval)
}
