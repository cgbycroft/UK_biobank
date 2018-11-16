

source("/well/ukbiobank/qcoutput.V2_QCed.sample-QC-testing/QC-Scripts/R/scripts/ukbbcolors.R")


## Routine are modified from SNPolisher
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
			 col=character(),cex=0.8,pch=c(4,19),
			 xrange=numeric(),yrange=numeric(),main.text="",
			 add=FALSE,add.ellipse=TRUE,axes=TRUE,
			 orig.scale=FALSE,axis.col=NULL,
			 xlab.text="Contrast:  log2(A/B)",
			 ylab.text="Strength:  log2(A*B)/2") {
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
    ylba.text = "B intensity"
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
    col.ellipse0 = geno.cols.lighter[2]
    col.ellipse1 = geno.cols.lighter[3]
    col.ellipse2 = geno.cols.lighter[4]
    lwd.ellipse = 3
  } else {
    color.indivs = col
    col.ellipse0 = "black"
    col.ellipse1 = "black"
    col.ellipse2 = "black"
    lwd.ellipse = 2
  }
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
  main.text = paste(main.text,"\n#samples = ",n," (#missing = ",m,")",sep="")
  ## sort.list sorts alphabetically by group name (i.e. color name)
  ##order = sort.list(color.indivs,decreasing = TRUE)
  ## but I think that the following will work to sort by number of occurrences
  order = sort.list(table(color.indivs)[color.indivs],decreasing = TRUE)
  if (!add) {
    plot(logratio[order],strength[order],col=color.indivs[order],
	 pch=char.indivs[order],cex=cex,bty="n",main=main.text,lwd=3,
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
    points(logratio[order],strength[order],col=color.indivs[order],
	   pch=char.indivs[order],cex=cex,lwd=3)
  }
  return(data.frame(xrange=xrange,yrange=yrange))
}

## Routines to read calls/intensities/snp-posteriors from packed binary files
unpack.probeset.data <- function(datadir,pid,is.autosomal,pname=NULL) {
  if (is.logical(is.autosomal)) {
    if (is.autosomal) { ProbeSet = "autosome" }
    else              { ProbeSet = "sexchrom" }
  } else {
    stop("unpack.probeset.data(datadir,pid,is.autosomal)")
  }
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
  seek(Binary,offset,origin="current")
  id = readBin(Binary,character(),n=1)
  id = sub(" *","",id)
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
