
myManhattan <- function (x, chr = "CHR", bp = "BP", p = "P", snp = "SNP", col = "black",chromCols=c("black","darkblue"), ymax = NULL,xlimits=NULL,suggestiveline = -log10(1e-05), genomewideline = -log10(5e-08), logtransform=TRUE,
    highlight = NULL,highlightCols = NULL,highlightPCH=20,highlightCEX=0.8,axisMB=TRUE,printAxes=TRUE,...)
{
    P = index = NULL
    if (!(chr %in% names(x))) 
        stop(paste("Column", chr, "not found!"))
    if (!(bp %in% names(x))) 
        stop(paste("Column", bp, "not found!"))
    if (!(p %in% names(x))) 
        stop(paste("Column", p, "not found!"))
    if (!(snp %in% names(x))) 
        warning(paste("No SNP column found. OK unless you're trying to highlight."))
    if (!is.numeric(x[[chr]])) 
        stop(paste(chr, "column should be numeric. Do you have 'X', 'Y', 'MT', etc? If so change to numbers and try again."))
    if (!is.numeric(x[[bp]])) 
        stop(paste(bp, "column should be numeric."))
    if (!is.numeric(x[[p]])) 
        stop(paste(p, "column should be numeric."))
    d = data.frame(CHR = x[[chr]], BP = x[[bp]], P = x[[p]])
    d$col = col
    if (!is.null(x[[snp]])) 
        d = transform(d, SNP = x[[snp]])    
    if(logtransform) d = subset(d[order(d$CHR, d$BP), ], (P >= 0 & P <= 1 & !is.na(P)))
    d$logp = -log10(d$P)
    if(!logtransform) d$logp = d$P
    d$pos = NA
    d$shape = 20
    if (!is.null(ymax)) {
        print( paste0("Capping ",sum(d$logp > ymax)," values above ",ymax))
        d$shape[d$logp > ymax] = 17        
        d$logp[d$logp > ymax] = ymax
    }
    if (is.null(ymax)) {
        ymax = ceiling(max(d$logp))
        message("Ymax will be set automatically based on most significant SNP")
    }
    else if (!is.numeric(ymax)) {
        ymax = ceiling(max(d$logp))
        warning("non-numeric ymax argument.")
    }
    else if (ymax < 0) {
        ymax = ceiling(max(d$logp))
        warning("negative ymax argument.")
    }
    else if (is.null(ymax)) {
        ymax = ceiling(max(d$logp))
        message(paste("Using", ymax, "as ymax"))
    }
    message(paste("Using", ymax, "as ymax"))
    d$index = NA
    ind = 0
    for (i in unique(d$CHR)) {
        ind = ind + 1
        d[d$CHR == i, ]$index = ind
    }
    nchr = length(unique(d$CHR))
    if (nchr == 1) {
        
        if(axisMB) {
            xlabel = paste("Chromosome", unique(d$CHR), "position (Mb)")
            d$BP = d$BP/1000000
        } else {
            xlabel = paste("Chromosome", unique(d$CHR), "position")
        }
            d$pos = d$BP
        
        #print(range(d$pos))      
        ticks = floor(length(d$pos))/2 + 1
        labs = ticks
    }
    else {
        lastbase = 0
        ticks = NULL
        for (i in unique(d$index)) {
            if (i == 1) {
                d[d$index == i, ]$pos = d[d$index == i, ]$BP
            }
            else {
                lastbase = lastbase + tail(subset(d, index == 
                  i - 1)$BP, 1)
                d[d$index == i, ]$pos = d[d$index == i, ]$BP + 
                  lastbase
            }
            ticks = c(ticks, d[d$index == i, ]$pos[floor(length(d[d$index == 
                i, ]$pos)/2) + 1])
        }
        xlabel = "Chromosome"
        labs = unique(d$CHR)
    }
        
    #xlimits=c(floor(max(d$pos) * -0.03),ceiling(max(d$pos) * 1.03))
    if (is.null(xlimits)){
      print("small region")
       xlimits <- c(floor(max(d$pos) * -0.03),ceiling(max(d$pos) * 1.03))
      #xlimits=c(min(d$pos)-floor(diff(range(d$pos))*0.03),ceiling(max(d$pos) * 1.03))      
      }
    
    #print(xlimits)
    ylabel = expression(-log[10](italic(p)))
    if(!printAxes) {
     xlabel = ylabel = NA   
    }
    ymin = min(0,min(d$logp) )
    plot(NULL, yaxt = "n", xaxt = "n", bty = "n", xaxs = "i", yaxs = "i", 
        xlim = xlimits, ylim = c(ymin, ymax + ymax/100), xlab = xlabel, 
        ylab = ylabel, las = 1, pch = 20, 
         ...)
    if(printAxes){
        
        axis(2, ...)
        
        if (nchr == 1) {
            axis(1, ...)
        }
        else {
            axis(1, at = ticks, labels = labs, ...)
        }
    }
    #if(is.null(col)) col = rep(c("gray10","gray60"), max(d$CHR))

    # print horizontal lines
    if (suggestiveline) 
        abline(h = suggestiveline, col = "blue",lwd=2, xpd=FALSE)
    if (genomewideline) 
        abline(h = genomewideline, col = "red",lwd=2, xpd=FALSE)

    # print p-value points
    if (nchr == 1) {
        with(d, points(pos, logp, pch = shape,col = col,xpd=NA,...))
    }
    else {
        icol = 1
        for (i in unique(d$index)) {
            with(d[d$index == unique(d$index)[i], ], points(pos, 
                logp, col = chromCols[icol%%2+1], pch = shape,xpd=NA,...))
            icol = icol + 1
        }
    }
    
    # print any highlighted points (on top of existing ones...)
    if (!is.null(highlight)) {
        if (any(!(highlight %in% d$SNP))) 
            warning("You're trying to highlight SNPs that don't exist in your results. These will not be printed.")
        d.highlight = d[which(d$SNP %in% highlight), ]
        
        if (!is.null(highlightCols)) highlightColours = highlightCols[match(d.highlight$SNP,highlight)] else highlightColours = "green3"

        with(d.highlight, points(pos, logp, col = highlightColours, pch = highlightPCH,xpd=NA,
            ...))
    }
}
