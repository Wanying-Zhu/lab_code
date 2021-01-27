#qq.chisq Quantile-quantile plot for chi-squared tests
#Description
#This function plots ranked observed chi-squared test statistics against the corresponding expected
#order statistics. It also estimates an inflation (or deflation) factor, lambda, by the ratio of the trimmed
#means of observed and expected values. This is useful for inspecting the results of whole-genome
#association studies for overdispersion due to population substructure and other sources of bias or
#confounding.
######## Added by Anna P. 
#  In addition to chi-squared test statistic, this function can plot
#  log-transformed P-values from any test, not necessarily a chi-squared
#  one. It uses the fact that -2*log(P-value), the NATURAL log, has a
#  chi-squared distribution with 2 df. The program determines whether the
#  input data is a list of P-values or a list of chi-squared statistic
#  values on the basis of a keyword in the axis or the subtitle labels. So
#  either one of them MUST contain the character sequence "p-val" or "chi"
#  (case INsensitive). The resulting plot is transformed to match
#  log10(P-value) 
######## The end of additional feature
#Usage
#qq.chisq(x, df=1, x.max, main="QQ plot",
#sub=paste("Expected distribution: chi-squared (",df," df)", sep=""),
#xlab="Expected", ylab="Observed",
#conc=c(0.025, 0.975), overdisp=FALSE, trim=0.5,
#slope.one=FALSE, slope.lambda=FALSE,
#thin=c(0.25,50), oor.pch=24, col.shade="gray", ...)
#Arguments
#x A vector of observed chi-squared test values
#df The degreees of freedom for the tests
#x.max If present, truncate the observed value (Y) axis here
#main The main heading
#sub The subheading
#xlab x-axis label (default "Expected")
#ylab y-axis label (default "Observed")
#conc Lower and upper probability bounds for concentration band for the plot. Set this
#to NA to suppress this
#overdisp If TRUE, an overdispersion factor, lambda, will be estimated and used in calculating
#concentration band
#trim Quantile point for trimmed mean calculations for estimation of lambda. Default
#is to trim at the median
#slope.one Is a line of slope one to be superimpsed?
#slope.lambda Is a line of slope lambda to be superimposed?
#thin A pair of numbers indicating how points will be thinned before plotting (see
#Details). If NA, no thinning will be carried out
#oor.pch Observed values greater than x.max are plotted at x.max. This argument sets
#the plotting symbol to be used for out-of-range observations
#col.shade The colour with which the concentration band will be filled
#... Further graphical parameter settings to be passed to points()

#Details
#To reduce plotting time and the size of plot files, the smallest observed and expected points are
#thinned so that only a reduced number of (approximately equally spaced) points are plotted. The
#precise behaviour is controlled by the parameter thin, whose value should be a pair of numbers.
#The first number must lie between 0 and 1 and sets the proportion of the X axis over which thinning
#is to be applied. The second number should be an integer and sets the maximum number of points
#to be plotted in this section.
#The "concentration band" for the plot is shown in grey. This region is defined by upper and lower
#probability bounds for each order statistic. The default is to use the 2.5 Note that this is not a
#simultaneous confidence region; the probability that the plot will stray outside the band at some
#point exceeds 95
#When required, the dispersion factor is estimated by the ratio of the observed trimmed mean to its
#expected value under the chi-squared assumption.
#Value
#The function returns the number of tests, the number of values omitted from the plot (greater than
#x.max), and the estimated dispersion factor, lambda.
#Note
#All tests must have the same number of degrees of freedom. If this is not the case, I suggest
#transforming to p-values and then plotting -2log(p) as chi-squared on 2 df.
#Author(s)
#David Clayton hdavid.clayton@cimr.cam.ac.uk
#References
#Devlin, B. and Roeder, K. (1999) Genomic control for association studies. Biometrics, 55:997-1004
# from snpMatrix
# http://www.bioconductor.org/packages/bioc/html/snpMatrix.html
 qq.chisq <-
  function(x, df=1, x.max,
    main="QQ plot",plotType="pval",
    sub=paste("Expected distribution: chi-squared (",df," df)", sep=""),
    xlab="Expected", ylab="Observed",
    conc=c(0.025, 0.975), dot.col="black", overdisp=FALSE, trim=0.5, 
    slope.one=T, slope.lambda=FALSE,
    thin=c(0.25,100), oor.pch=24, col.shade="gray", ofname="qqchi.pdf",
    h=6,w=6,printpdf=F,...) {

# Added by Anna P. to allow for plotting log(P-value) instead of 
# 2*log(P-value) to make it compatible with chisq distribution with 2 df
    if (plotType == "pval") {
      x <- -2*log(x)
      scale <- 0.5/log(10)
      x.max <- x.max/scale
      df <- 2
    } else {
      scale <- 1
    }

    # Function to shade concentration band
   
    shade <- function(x1, y1, x2, y2, color=col.shade) {
      n <- length(x2)
      polygon(c(x1, x2[n:1]), c(y1, y2[n:1]), border=NA, col=color)
    }
   
    # Sort values and see how many out of range
   
    obsvd <- sort(x, na.last=NA)
    N <- length(obsvd)
    if (missing(x.max)) {
      Np <- N
    }
    else {
      Np <- sum(obsvd<=x.max)
    }
    if(Np==0)
      stop("Nothing to plot")

    # Expected values
   
    if (df==2) {
      expctd <- 2*cumsum(1/(N:1))
    }
    else {
      expctd <- qchisq(p=(1:N)/(N+1), df=df)
    }

    # Concentration bands
   
    if (!is.null(conc)) {
      if(conc[1]>0) {
        e.low <- qchisq(p=qbeta(conc[1], 1:N, N:1), df=df)
      }
      else {
        e.low <- rep(0, N)
      }
      if (conc[2]<1) {
        e.high <- qchisq(p=qbeta(conc[2], 1:N, N:1), df=df)
      }
      else {
        e.high <- 1.1*rep(max(x),N)
      }
    }
   
    # Plot outline
   
    if (Np < N)
      top <- x.max
    else
      top <- obsvd[N]
    right <- expctd[N]
    if (printpdf) {pdf(ofname,h,w)}
    plot(c(0, right*scale), c(0, top*scale), type="n", xlab=xlab, ylab=ylab, bty="l",
         main=main, sub=sub, cex.lab=1.25,
	 cex.axis=1.1,ylim=c(0,x.max*scale), xlim=c(0,9))
   
    # Thinning
   
    if (is.na(thin[1])) {
      show <- 1:Np
    }
    else if (length(thin)!=2 || thin[1]<0 || thin[1]>1 || thin[2]<1) {
      warning("invalid thin parameter; no thinning carried out")
      show <- 1:Np
    }
    else {
      space <- right*thin[1]/floor(thin[2])
      iat <- round((N+1)*pchisq(q=(1:floor(thin[2]))*space, df=df))
      if (max(iat)>thin[2])
        show <- unique(c(iat, (1+max(iat)):Np))
      else
        show <- 1:Np
    }
    Nu <- floor(trim*N)
#    Nl <- floor(0.3*N)
    if (Nu>0)
      lambda <- mean(obsvd[1:Nu])/mean(expctd[1:Nu])
#      lambda <- mean(obsvd[Nl:Nu])/mean(expctd[Nl:Nu])
    if (!is.null(conc)) {
      if (Np<N)
        vert <- c(show, (Np+1):N)
      else
        vert <- show
      if (overdisp)
        shade(expctd[vert]*scale, lambda*e.low[vert]*scale,
              expctd[vert]*scale, lambda*e.high[vert]*scale)
      else
        shade(expctd[vert]*scale, e.low[vert]*scale, expctd[vert]*scale, e.high[vert]*scale)
    }
    points(expctd[show]*scale, obsvd[show]*scale,cex=0.5,col=dot.col, ...)
    # Overflow
    if (Np<N) {
      over <- (Np+1):N
      points(expctd[over]*scale, rep(x.max, N-Np)*scale, pch=oor.pch,col=dot.col)
    }
    # Lines
    line.types <- c("solid", "dashed", "dotted")
    key <- NULL
    txt <- NULL
    if (slope.one) {
      key <- c(key, line.types[1])
      txt <- c(txt, "y = x")
      abline(a=0, b=1, lty=line.types[1])
    }
    if (slope.lambda && Nu>0) {
      key <- c(key, line.types[2])       
      txt <- c(txt, paste("y = ", format(lambda, digits=4), "x", sep=""))
      if (!is.null(conc)) {
        if (Np<N)
          vert <- c(show, (Np+1):N)
        else
          vert <- show
      }   
      abline(a=0, b=lambda, lty=line.types[2])
    }
    if (printpdf) {dev.off()}
    # Returned value
   
#    if (!is.null(key))
#       legend(0, top*scale, legend=txt, lty=key, bty="n")
    c(N=N, omitted=N-Np, lambda=lambda)

  }
