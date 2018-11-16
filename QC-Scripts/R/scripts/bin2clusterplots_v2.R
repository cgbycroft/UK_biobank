
source("/well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/R/scripts/ukbbcolors.R")


## Routine are modified from SNPolisher.
## If you want exactly 85% of the probability density to be inside the ellipse, then set c=sqrt(qchisq(0.85,df = 2)) = 1.947881

ellipse.info <- function(x0,vx,y0,vy,cov,c=2) {
  theta = 0.5 * atan2(2*cov, vx-vy)
  sint = sin(theta)
  cost = cos(theta)
  np = 100
  a = c*sqrt(vx*cost*cost + vy*sint*sint + cov*2*sint*cost)
  b = c*sqrt(vx*sint*sint + vy*cost*cost - cov*2*sint*cost)
  alpha = 2 * pi * (0:np)/np
  sint = sin(theta)
  cost = cos(theta)
  sina = sin(alpha)
  cosa = cos(alpha)
  xpts = x0 + a * cosa * cost - b * sina * sint
  ypts = y0 + a * cosa * sint + b * sina * cost
  return(data.frame(xpts,ypts))
}

cluster.info <- function(params,c) {
  df = NULL
  if (is.numeric(params)) {
    m = length(params)/3
    ## Why 7? This is how the parameters of a posterior cluster are stored in snp-posteriors.bin
    if (m!=7) { stop("m!=7") }
    params = params[c*m+1:m]
    ## I will write the information in slightly different format: c(mu,sigma,n)                           
    ## where mu = c(x,y), sigma = c(sigma11,sigma12,sigma21,sigma22) and n is NObsMean
    df = data.frame(x0=params[1],vx=params[3],
      y0=params[2],vy=params[6],
      cov=params[4],n=params[7])
  }
  return(df)
}

cluster.plot <- function(data,
			 exclude.indivs=numeric(),
			 col=character(),cex=0.6,pch=c(4,19),
			 xrange=numeric(),yrange=numeric(),main.text="",
			 add=FALSE,add.ellipse=TRUE,axes=TRUE,
			 orig.scale=FALSE,axis.col=NULL,
			 xlab.text="Contrast:  log2(A/B)",
			 ylab.text="Strength:  log2(A*B)/2",Alpha=1,showCounts=TRUE,showMaf=FALSE) {
    posteriors = NULL
    posteriors1 = NULL
    if (!is.null(data$intensities)) {
        intensities = data$intensities
        n = length(intensities)/2
        logratio = intensities[2*(1:n)-1]
        strength = intensities[2*(1:n)]    
    } else if (!is.null(data$logratio) && !is.null(data$strength)) {
        logratio = data$logratio
        strength = data$strength
        n = length(logratio)
    } else {
        stop("cluster.plot : Provide either intensities or logratio+strength")
    }
    if (!is.null(data$posteriors)) { posteriors = data$posteriors }
    if (!is.null(data$posteriors1)) { posteriors1 = data$posteriors1 }
    if (!is.null(data$calls)) {
        calls = data$calls
        calls[is.na(calls)] = -1
    } else {
        stop("cluster.plot : Provide calls")
    }
    if (orig.scale) {
        a = 2^(strength + logratio/2)
        b = 2^(strength - logratio/2)
        logratio = a
        strength = b
        xlab.text = "A intensity"
        ylab.text = "B intensity"
        add.ellipse = FALSE
    }
    cluster0 = NULL; cluster1 = NULL; cluster2 = NULL;
    ellipse0 = NULL; ellipse1 = NULL; ellipse2 = NULL;
    ellipse0.1 = NULL; ellipse1.1 = NULL; ellipse2.1 = NULL;
    cluster0.1 = NULL; cluster1.1 = NULL; cluster2.1 = NULL;
    if (!is.null(posteriors)) {
        cluster0 = cluster.info(posteriors,c=0)
        cluster1 = cluster.info(posteriors,c=1)
        cluster2 = cluster.info(posteriors,c=2)
        if ((cluster0$vx>0)&(cluster1$vx)&(cluster2$vx)&
            (cluster0$vy>0)&(cluster1$vy)&(cluster2$vy)) {
            ellipse0 = ellipse.info(cluster0$x0,cluster0$vx,cluster0$y0,cluster0$vy,cluster0$cov)
            ellipse1 = ellipse.info(cluster1$x0,cluster1$vx,cluster1$y0,cluster1$vy,cluster1$cov)
            ellipse2 = ellipse.info(cluster2$x0,cluster2$vx,cluster2$y0,cluster2$vy,cluster2$cov)
        }
    }
    if (!is.null(posteriors1)) {
        cluster0.1 = cluster.info(posteriors1,c=0)
        cluster1.1 = cluster.info(posteriors1,c=1)
        cluster2.1 = cluster.info(posteriors1,c=2)
        if ((cluster0.1$vx>0)&(cluster1.1$vx)&(cluster2.1$vx)&
            (cluster0.1$vy>0)&(cluster1.1$vy)&(cluster2.1$vy)) {
            ellipse0.1 = ellipse.info(cluster0.1$x0,cluster0.1$vx,cluster0.1$y0,cluster0.1$vy,cluster0.1$cov)
            ellipse1.1 = ellipse.info(cluster1.1$x0,cluster1.1$vx,cluster1.1$y0,cluster1.1$vy,cluster1.1$cov)
            ellipse2.1 = ellipse.info(cluster2.1$x0,cluster2.1$vx,cluster2.1$y0,cluster2.1$vy,cluster2.1$cov)
        }
    }
    if (!length(col)) {
        color.indivs = geno.cols[calls+2]
        #col.ellipse0 = geno.cols.lighter[2]
        #col.ellipse1 = geno.cols.lighter[3]
        #col.ellipse2 = geno.cols.lighter[4]
        #lwd.ellipse = 3
        col.ellipse0 = col.ellipse1 = col.ellipse2 = "black"
        lwd.ellipse = 1.2
        
    } else {
        color.indivs = col
        col.ellipse0 = "black"
        col.ellipse1 = "black"
        col.ellipse2 = "black"
        lwd.ellipse = 1.2
    }

    if(sum(calls!=(-1))>0) color.indivs[calls!=(-1)] = add.alpha(color.indivs[calls!=(-1)],Alpha)
    
    if (length(pch)==n) {
        char.indivs = pch
    } else if (length(pch)==1) {
        char.indivs = rep(pch,n)
    } else if (length(pch)==2) {
        char.indivs = pch[1+1*(calls==0)+1*(calls==1)+1*(calls==2)]
    } else if (length(pch)==4) {
        char.indivs = pch[1+1*(calls==0)+2*(calls==1)+3*(calls==2)]
    } else {
        char.indivs = rep(1,n)
    }
    if (length(exclude.indivs)) {
        calls = calls[-exclude.indivs]
        logratio = logratio[-exclude.indivs]
        strength = strength[-exclude.indivs]
        char.indivs = char.indivs[-exclude.indivs]
        color.indivs = color.indivs[-exclude.indivs]
    }
    
    
    if (!length(xrange)) {
        xrange = range(c(logratio,ellipse0$xpts,ellipse1$xpts,ellipse2$xpts,
            ellipse0.1$xpts,ellipse1.1$xpts,ellipse2.1$xpts))
    }
    if (!length(yrange)) {
        yrange = range(c(strength,ellipse0$ypts,ellipse1$ypts,ellipse2$ypts,
            ellipse0.1$ypts,ellipse1.1$ypts,ellipse2.1$ypts))
    }
    n = length(calls)
    m = sum(calls==-1)
    if (is.null(axis.col)) { axis.col = "black" }
    
    
    if (!length(logratio)) {
        if (!add) {
            plot(0,0,type="n",bty="n",main="",axes=FALSE,
                 col.axis=axis.col,col.lab=axis.col,
                 col.main=axis.col,col.sub=axis.col,
                 xlab=xlab.text,xlim=xrange,
                 ylab=ylab.text,ylim=yrange)
            if (axes) {
                axis(1,col=axis.col,col.axis=axis.col,col.lab=axis.col,col.ticks=axis.col)
                axis(2,col=axis.col,col.axis=axis.col,col.lab=axis.col,col.ticks=axis.col)
            } else {
                box( )
            }
        }
        return(data.frame(xrange=xrange,yrange=yrange))
    }

    callsCounts = paste(c(sum(calls==0),sum(calls==1),sum(calls==2),m),collapse="|")
    if(sum(calls%in%c(0,1,2))>0){
        maf = (2*sum(calls==0)+sum(calls==1))/(2*sum(calls%in%c(0,1,2)))
        mac = 2*sum(calls==0)+sum(calls==1)
#        print(table(calls))
        if(maf>0.5) {
            mac = 2*sum(calls==2)+sum(calls==1)
            maf = 1-maf
        }
    } else {
        # if there are no calls non-missing calls!
        maf=NA
        mac=NA
    }
    #main.text = paste(main.text,"\n#samples = ",n," (#missing = ",m,")",sep="")
    if(showCounts) main.text2 = paste0("n=",n,"; AA|AB|BB|nocall = ",callsCounts) else main.text2=paste0("#samples = ",n," (#no calls = ",m,")")
    if(showMaf%in%c(1,2)) {
        main.text2 = paste0("n = ",n,"; #miss = ",m,"; maf = ",round(maf,3))
        if(( !is.na(maf)&(maf<0.001) ) | (showMaf==2) ) main.text2 = paste0("n = ",n,"; #miss = ",m,"; mac = ",mac)
        
    } else {
        main.text2=paste0("#samples = ",n," (#no calls = ",m,")")
    }
    
    ## sort.list sorts alphabetically by group name (i.e. color name)
    ##order = sort.list(color.indivs,decreasing = TRUE)
    ## but I think that the following will work to sort by number of occurrences
    order = sort.list(table(color.indivs)[color.indivs],decreasing = TRUE)

    if (!add) {
        plot(logratio[order],strength[order],col=color.indivs[order],bg=color.indivs[order],
             pch=char.indivs[order],cex=cex,bty="n",main=main.text,lwd=1.5,
             col.axis=axis.col,col.lab=axis.col,axes=FALSE,
             col.main=axis.col,col.sub=axis.col,
             xlab=xlab.text,xlim=xrange,
             ylab=ylab.text,ylim=yrange)
        if (axes) {
            axis(1,col=axis.col,col.axis=axis.col,col.lab=axis.col,col.ticks=axis.col)
            axis(2,col=axis.col,col.axis=axis.col,col.lab=axis.col,col.ticks=axis.col)
        } else {
            box( )
        }
        if (add.ellipse) {
            if (!is.null(ellipse0)&!is.null(ellipse1)&!is.null(ellipse2)) {
                lines(ellipse0$xpts,ellipse0$ypts,type="l",lwd=lwd.ellipse,col=col.ellipse0)
                lines(ellipse1$xpts,ellipse1$ypts,type="l",lwd=lwd.ellipse,col=col.ellipse1)
                lines(ellipse2$xpts,ellipse2$ypts,type="l",lwd=lwd.ellipse,col=col.ellipse2)
            }
            if (!is.null(ellipse0.1)&!is.null(ellipse1.1)&!is.null(ellipse2.1)) {
                lines(ellipse0.1$xpts,ellipse0.1$ypts,type="l",lwd=lwd.ellipse,col=col.ellipse0,lty=3)
                lines(ellipse1.1$xpts,ellipse1.1$ypts,type="l",lwd=lwd.ellipse,col=col.ellipse1,lty=3)
                lines(ellipse2.1$xpts,ellipse2.1$ypts,type="l",lwd=lwd.ellipse,col=col.ellipse2,lty=3)
            }
        }
    } else {
        points(logratio[order],strength[order],col=color.indivs[order],bg=color.indivs[order],
               pch=char.indivs[order],cex=cex,lwd=1.5)
    }
    mtext(side=3,line=0.5,text=main.text2,cex=par()$cex.main/1.5)
#    mtext(side=3,line=0.5,text=main.text2,cex=par()$cex.main/3) #### <=== TEMP!

    return(data.frame(xrange=xrange,yrange=yrange))
}


#############
## Theses are for binary files created in qc pipeline (from affymetrix intensity and calls files).
## Routines to read calls/intensities/snp-posteriors from packed binary files
#############

## Routines to read calls/intensities/snp-posteriors from packed binary files
unpack.probeset.data <- function(datadir,pid,is.autosomal,pname=NULL) {

    writeLines("sheep")
    
  if (is.logical(is.autosomal)) {
    if (is.autosomal) { ProbeSet = "autosome" }
    else              { ProbeSet = "sexchrom" }
  } else {
    stop("unpack.probeset.data(datadir,pid,is.autosomal)")
  }

  writeLines(paste(datadir,'/AxiomGT1.',ProbeSet,'.calls.bin',sep=''))
  writeLines(paste(datadir,'/AxiomGT1.',ProbeSet,'.intensities.bin',sep=''))
  writeLines(paste(datadir,'/AxiomGT1.',ProbeSet,'.snp-posteriors.bin',sep=''))
  
  calls = unpack.calls(paste(datadir,'/AxiomGT1.',ProbeSet,'.calls.bin',sep=''),pid,pname=pname)
  intensities = unpack.intensities(paste(datadir,'/AxiomGT1.',ProbeSet,'.intensities.bin',sep=''),pid,pname=pname)
  posteriors = unpack.snp.posterior(paste(datadir,'/AxiomGT1.',ProbeSet,'.snp-posteriors.bin',sep=''),pid,pname=pname)
  if (is.null(calls) || is.null(intensities) || is.null(posteriors)) {
    return(NULL)
  }
  if (ProbeSet=="sexchrom") {
    ## First diploid posterior, then haploid posterior
    posteriors1 = posteriors[22:42]
    posteriors = posteriors[1:21]
  } else {
    posteriors1 = NULL
  }
  data = list(calls=calls,intensities=intensities,posteriors=posteriors,posteriors1=posteriors1)
  return(data)
}

unpack.calls <- function(Datafile,pid,pname=NULL) {
    
    Binary = file(Datafile,"rb")
    
    values = readBin(Binary,numeric(),n=3)
    RowNum = as.integer(values[1])
    RecNum = as.integer(values[2])
    RecSize = as.integer(values[3])
    
    offset = (pid-1)*(RecSize*RecNum+16)

    writeLines(paste(RowNum,RecNum,RecSize,"offset:",offset,"pid:",pid,"pname:",pname))

    seek(Binary,offset,origin="current")
    id = readBin(Binary,character(),n=1)
    id = sub(" *","",id)

    writeLines(paste("id:",id))
    
    calls = readBin(Binary,what="int",size=1,n=RecNum)
    close(Binary)
    if (!is.null(pname) && (pname!=id)) {
        calls = NULL
    }
    
    return(calls)
}

unpack.snp.posterior <- function(Datafile,pid,pname=NULL) {
  Binary = file(Datafile,"rb")
  values = readBin(Binary,numeric(),n=3)
  RowNum = as.integer(values[1])
  RecNum = as.integer(values[2])
  RecSize = as.integer(values[3])
  offset = (pid-1)*(RecSize*RecNum+16)
  tempi = seek(Binary,offset,origin="current")
  id = readBin(Binary,character(),n=1)
  id = sub(" *","",id)
  posterior = readBin(Binary,numeric(),size=RecSize,n=RecNum)
  close(Binary)
  
  if (!is.null(pname) && (pname!=id)) {
    posterior = NULL
  }

  return(posterior)
}

unpack.intensities <- function(Datafile,pid,pname=NULL) {
    Binary = file(Datafile,"rb")
    values = readBin(Binary,numeric(),n=3)
    RowNum = as.integer(values[1])
    RecNum = as.integer(values[2])
    RecSize = as.integer(values[3])
    offset = (pid-1)*(2*RecSize*RecNum+16)
    seek(Binary,offset,origin="current")
    id = readBin(Binary,character(),n=1)
    id = sub(" *","",id)
    intensities = readBin(Binary,numeric(),size=RecSize,n=2*RecNum)
    close(Binary)
    if (!is.null(pname) && (pname!=id)) {
        intensities = NULL
    }

    return(intensities)
}



#############
## Theses are for binary files created for data release output.
## Routines to read calls/intensities/snp-posteriors from packed binary files
#############
                                        # ../expt/V2_QCed.export/src/export_intensity/probe_binary_intensity.pl 


# where are the calls/intensity/posterior files?
################
callsFiles = "/well/ukbiobank/expt/V2_QCed.export/data/calls_export/v3/V2_QCed.export.UKBiLEVEAX_b1-b11.Batch_b001_b095"
famFile=NULL
intensityFiles = "/well/ukbiobank/expt/V2_QCed.export/data/export_intensity/v2/V2_QCed.export.intensity"
posteriorFiles = "/well/ukbiobank/expt/V2_QCed.export/data/calls_export/v3/V2_QCed.export.UKBiLEVEAX_b1-b11.Batch_b001_b095"
################


getPID <- function(bimFiles=callsFiles,SNPs){
    # find line number of a set of snps in bim/intensity/calls files
    pidsRaw = system( paste0("grep -E -wn \'",paste(SNPs,collapse="|"),"\' ",bimFiles,"*.bim | awk '{print $1,$2}'") , intern=T)
    pids = t(sapply(pidsRaw,function(x) str_split(x,":| ")[[1]],USE.NAMES=FALSE))

    headerStuff = c("bim","line","chrom","pname")
    if(length(pids)==0) pids = as.data.frame(matrix(vector(),0,length(headerStuff) ))    
    colnames(pids) = headerStuff
    return(pids)
}

unpack.calls2 <- function(Datafile,famFile=NULL,pid,readFam=FALSE) {

    # NOTE: pid must be the line number of the desired snp in the order of the *.bim file.

    if(is.null(famFile) ) fam <- gsub("\\.bed",".fam",Datafile) else fam = famFile
    
    ## calculate number of samples in file
    ni <- as.numeric( system(paste0("wc -l ",fam," | cut -f1 -d' '"),intern=TRUE) )

    ## open bed file and check its magic number    
    bed <- file(Datafile, open = "rb")      
    magic <- readBin(bed, "raw", 3)
    
    if (!all(magic[1] == "6c", magic[2] == "1b", magic[3] == "01"))
        stop("Wrong magic number for bed file; should be -- 0x6c 0x1b 0x01 --.")


    ## block size in bytes: (number of individuals)/4, to nearest byte
    bsz <- ceiling(ni/4)
    
                                        # calculate offset for this variant; go there
    v <- pid-1
    offset <- v*bsz+3
    seek(bed, offset, "start")
                                        # read the genotypes
    geno.raw <- as.logical(rawToBits(readBin(bed, "raw", bsz)))    
    j <- seq(1,length(geno.raw),2)
    
    ## express genotypes as allele dosage (0,1,2) where it counts the number of copies of the second allele in the corresponding .bim file. This should match the order of alleles in the intensity files.    
    geno <- geno.raw[ j ] + geno.raw[ j+1 ]
    
    ## recall that 0/1 is het, but 1/0 is missing
    geno[ geno.raw[ j ] == 1 & geno.raw[ j+1 ] == 0 ] <- NA
                                        # NOTE: exclude the last zeros which are extra padding based on N/4. See https://www.cog-genomics.org/plink/1.9/formats#bed
    calls = geno[1:ni]

    if(readFam) {
        f = read.table(fam,header=FALSE,stringsAsFactors=FALSE)
        names(geno) = f[,2]
    }
        
    close(bed)
    
    return(calls)
    
}


unpack.snp.posterior2 <- function(datadir,pid,pname=NULL,is.autosomal=TRUE) {

    
    if (is.logical(is.autosomal)) {
        if (is.autosomal) { ProbeSet = "autosome" }
        else              { ProbeSet = "sexchrom" }
    } else {
        stop("unpack.probeset.data(datadir,pid,is.autosomal)")
    }

    Datafile=paste0(datadir,'/AxiomGT1.',ProbeSet,'.snp-posteriors.bin')

                                        # CURRENTLY THE SAME AS BEFORE!
    Binary = file(Datafile,"rb")
    values = readBin(Binary,numeric(),n=3)
    RowNum = as.integer(values[1])
    RecNum = as.integer(values[2])
    RecSize = as.integer(values[3])
    offset = (pid-1)*(RecSize*RecNum+16)
    tempi = seek(Binary,offset,origin="current")
    id = readBin(Binary,character(),n=1)
    id = sub(" *","",id)
    posteriors = readBin(Binary,numeric(),size=RecSize,n=RecNum)
    close(Binary)
    
    if (!is.null(pname) && (pname!=id)) {
        posteriors = NULL
    }

    
    if (ProbeSet=="sexchrom") {
        ## First diploid posterior, then haploid posterior
        posteriors1 = posteriors[22:42]
        posteriors = posteriors[1:21]
    } else {
        posteriors1 = NULL
    }
    
    return(list(posteriors=posteriors,posteriors1=posteriors1))
}


unpack.intensities2 <- function(Datafile,famFile=NULL,pid,readFam=FALSE) {

    if(is.null(famFile) ) fam <- gsub("\\.bin",".fam",Datafile) else fam = famFile

    ni <- as.numeric( system(paste0("wc -l ",fam," | cut -f1 -d' '"),intern=TRUE) )

    int <- file(Datafile, open = "rb")      

    bytes_per_snp=ni*4*2;
    num_fields=ni*2;
    
    bin_line=0
    int_byte_length=8; # length of bytes per sample
    offset=0;

    bin_seek = bytes_per_snp*( pid - 1 )
    seek(int, bin_seek,"start") # go to start of the right line.
    intRaw = readBin(int, numeric(), n=2*ni,size=4)
                                        # -1 values =NA
    intRaw[intRaw==-1] = NA
    
    close(int)

    
    # these are untransformed A and B intensities. Transform them!

    Aint = intRaw[2*(1:ni)-1]
    Bint = intRaw[2*(1:ni)]
    
    logratio = log2(Aint/Bint) # log2(A/B)
    strength = log2(Aint*Bint)/2   # log2(A*B)/2

    intensities = intRaw;
    intensities[2*(1:ni)-1] = logratio
    intensities[2*(1:ni)] = strength
    
    if(readFam) {
        f = read.table(fam,header=FALSE,stringsAsFactors=FALSE)
        return(list(intensities=intensities,samples=f))
    }
    
    return(intensities)    
}
