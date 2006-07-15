"abvnonpar"<- 
function(x = 0.5, data, k = NULL, epmar = FALSE, nsloc1 = NULL, nsloc2 = NULL,
         method = c("cfg","pickands","deheuvels","halltajvidi","tdo"),
         convex = FALSE, rev = FALSE, madj = 0, kmar = NULL, plot = FALSE,
         add = FALSE, lty = 1, lwd = 1, col = 1, blty = 3, blwd = 1,
         xlim = c(0,1), ylim = c(0.5,1), xlab = "t", ylab = "A(t)", ...)
{
    if(!is.null(k)) {
      return(abvnonparpot(x = x, data = data, k = k, epmar = epmar,
        convex = convex, rev = rev, kmar = kmar, plot = plot, add =
        add, lty = lty, lwd = lwd, col = col, blty = blty, blwd = blwd,
        xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab, ...))
    }
    if(mode(x) != "numeric" || any(x < 0, na.rm=TRUE) ||
       any(x > 1, na.rm=TRUE)) stop("invalid argument for `x'")
    if(epmar) {
      data <- apply(data, 2, rank, na.last = "keep")
      nasm <- apply(data, 2, function(x) sum(!is.na(x)))
      data <- data / rep(nasm+1, each = nrow(data))
      data <- -log(data)
    }
    else {
      if(is.null(kmar)) {
        if(!is.null(nsloc1)) {
            nsloc1 <- nsloc.transform(data, nsloc1)
            nslocmat1 <- cbind(1,as.matrix(nsloc1))
        }
        if(!is.null(nsloc2)) {
            nsloc2 <- nsloc.transform(data, nsloc2)
            nslocmat2 <- cbind(1,as.matrix(nsloc2))
        }
        # Transform to exponential margins
        mle.m1 <- frobgev(data[,1], nsloc = nsloc1)
        loc.mle.m1 <- mle.m1[grep("^loc", names(mle.m1))]
        if(is.null(nsloc1)) loc.mle.m1 <- rep(loc.mle.m1, nrow(data))
        else loc.mle.m1 <- nslocmat1 %*% loc.mle.m1
        mle.m1 <- cbind(loc.mle.m1, mle.m1["scale"], mle.m1["shape"])
        mle.m2 <- frobgev(data[,2], nsloc = nsloc2)
        loc.mle.m2 <- mle.m2[grep("^loc", names(mle.m2))]
        if(is.null(nsloc2)) loc.mle.m2 <- rep(loc.mle.m2, nrow(data))
        else loc.mle.m2 <- nslocmat2 %*% loc.mle.m2
        mle.m2 <- cbind(loc.mle.m2, mle.m2["scale"], mle.m2["shape"])
        data <- mtransform(data, list(mle.m1, mle.m2))       
        # End transform
      }
      else {
        if(!is.null(nsloc1) || !is.null(nsloc2))
            warning("ignoring `nsloc1' and `nsloc2' arguments")
        data <- mtransform(data, kmar)
      }
    }
    
    if(rev) data <- data[,2:1]
    data <- na.omit(data)
    if(plot || add) x <- seq(0, 1, length = 100)
    method <- match.arg(method)
    mpmin <- function(a,b) {
      a[a > b] <- b[a > b]
      a
    }
    mpmax <- function(a,b) {
      a[a < b] <- b[a < b]
      a
    }
    nn <- nrow(data)

    if((method == "pickands") && (madj < 0.5)) {
        if(!convex) {
            a <- numeric(length(x))
            for(i in 1:length(x))
                a[i] <- sum(mpmin(data[,1]/x[i], data[,2]/(1-x[i])))
            a <- nn / a
            a <- pmin(1, pmax(a, x, 1-x))
        }
        else {
            x2 <- seq(0, 1, length = 250)
            a <- numeric(250)
            for(i in 1:250)
                a[i] <- sum(mpmin(data[,1]/x2[i], data[,2]/(1-x2[i])))
            a <- nn / a
            a <- pmin(1, pmax(a, x2, 1-x2))
            inch <- chull(x2, a)
            a <- a[inch] ; x2 <- x2[inch]
            a <- approx(x2, a, xout = x, method="linear")$y
        }
    } 
    if((method == "deheuvels") ||
      ((method == "pickands") && (madj >= 0.5) && (madj < 1.5))) {
        if(!convex) {
            a <- numeric(length(x))
            for(i in 1:length(x))
                a[i] <- sum(mpmin(data[,1]/x[i], data[,2]/(1-x[i])))
            a <- nn / (a - x * sum(data[,1]) - (1-x) * sum(data[,2]) + nn)
            a <- pmin(1, pmax(a, x, 1-x))
        }
        else {
            x2 <- seq(0, 1, length = 250)
            a <- numeric(250)
            for(i in 1:250)
                a[i] <- sum(mpmin(data[,1]/x2[i], data[,2]/(1-x2[i]))) 
            a <- nn / (a - x2 * sum(data[,1]) - (1-x2) * sum(data[,2]) + nn)
            a <- pmin(1, pmax(a, x2, 1-x2))
            inch <- chull(x2, a)
            a <- a[inch] ; x2 <- x2[inch]
            a <- approx(x2, a, xout = x, method="linear")$y
        }
    }
    if(method == "cfg") {
        slm1 <- sum(log(data[,1]))
        slm2 <- sum(log(data[,2]))
        if(!convex) {
            a <- numeric(length(x))
            for(i in 1:length(x))
                a[i] <- sum(log(mpmax((1-x[i])*data[,1], x[i]*data[,2])))
            a <- (a - (1-x)*slm1 - x*slm2)/nn
            a <- pmin(1, pmax(exp(a), x, 1-x))
        }
        else {
            x2 <- seq(0, 1, length = 250)
            a <- numeric(250)
            for(i in 1:250)
                a[i] <- sum(log(mpmax((1-x2[i])*data[,1], x2[i]*data[,2])))
            a <- (a - (1-x2)*slm1 - x2*slm2)/nn
            a <- pmin(1, pmax(exp(a), x2, 1-x2))
            inch <- chull(x2, a)
            a <- a[inch] ; x2 <- x2[inch]
            a <- approx(x2, a, xout = x, method="linear")$y
        }
    }
    if(method == "tdo") {
        if(!convex) {
            a <- numeric(length(x))
            for(i in 1:length(x))
                a[i] <- sum(mpmin(x[i]/(1 + nn*data[,1]),
                                 (1-x[i])/(1 + nn*data[,2])))
            a <- 1 - a/(1 + log(nn))
            a <- pmin(1, pmax(a, x, 1-x))
        }
        else {
            x2 <- seq(0, 1, length = 250)
            a <- numeric(250)
            for(i in 1:250)
                a[i] <- sum(mpmin(x2[i]/(1 + nn*data[,1]),
                                 (1-x2[i])/(1 + nn*data[,2])))
            a <- 1 - a/(1 + log(nn))
            a <- pmin(1, pmax(a, x2, 1-x2))
            inch <- chull(x2, a)
            a <- a[inch] ; x2 <- x2[inch]
            a <- approx(x2, a, xout = x, method="linear")$y
        }
    } 
    if(method == "halltajvidi" ||
       ((method == "pickands") && (madj >= 1.5))) {
        sum1 <- sum(data[,1])
        sum2 <- sum(data[,2])
        if(!convex) {
            a <- numeric(length(x))
            for(i in 1:length(x))
                a[i] <- sum(mpmin(data[,1]/(sum1 * x[i]),
                                 data[,2]/(sum2 * (1 - x[i]))))
            a <- 1/a
            a <- pmin(1, pmax(a, x, 1-x))
        }
        else {
            x2 <- seq(0, 1, length = 250)
            a <- numeric(250)
            for(i in 1:250)
                a[i] <- sum(mpmin(data[,1]/(sum1 * x2[i]),
                                 data[,2]/(sum2 * (1 - x2[i]))))
            a <- 1 / a
            a <- pmin(1, pmax(a, x2, 1-x2))
            inch <- chull(x2, a)
            a <- a[inch] ; x2 <- x2[inch]
            a <- approx(x2, a, xout = x, method="linear")$y
        }
    } 
    if(plot || add) {
        bvdepfn(x, a, add, lty, lwd, col, blty, blwd, xlab, ylab, xlim, ylim, ...)  
        return(invisible(list(x = x, y = a)))
    } 
    a
}

"qcbvnonpar"<- 
function(p = seq(0.75, 0.95, 0.05), data, epmar = FALSE, nsloc1 = NULL,
         nsloc2 = NULL, mint = 1, method = c("cfg","pickands","deheuvels",
         "halltajvidi", "tdo"), convex = FALSE, madj = 0, kmar = NULL,
         plot = FALSE, add = FALSE, lty = 1, lwd = 1, col = 1, xlim =
         range(data[,1], na.rm = TRUE), ylim = range(data[,2], na.rm =
         TRUE), xlab = colnames(data)[1], ylab = colnames(data)[2], ...)
{
    if(mode(p) != "numeric" || p <= 0 || p >= 1)
      stop("`p' must be a vector of probabilities")
    nxv <- 100
    x <- seq(0, 1, length = nxv)
    ax <- abvnonpar(x = x, data = data, epmar = epmar, nsloc1 = nsloc1,
      nsloc2 = nsloc2, method = method, convex = convex, madj = madj,
      kmar = kmar, plot = FALSE)
    np <- length(p)
    qct <- list()
    p <- p^mint
    if(add) {
      xlim <- par("usr")[1:2]
      ylim <- par("usr")[1:2]
      if(par("xlog")) xlim <- 10^xlim
      if(par("ylog")) ylim <- 10^ylim
    }
    for(i in 1:np) {
      qct[[i]] <- -cbind(x/ax * log(p[i]), (1-x)/ax * log(p[i]))
      if(epmar) {
        qct[[i]] <- cbind(quantile(data[,1], probs = exp(-qct[[i]][,1]),
          na.rm = TRUE), quantile(data[,2], probs = exp(-qct[[i]][,2]),
          na.rm = TRUE))
      }
      else {
        if(is.null(kmar)) {
          # Transform from exponential margins
          mle.m1 <- frobgev(data[,1], nsloc = nsloc1)
          mle.m2 <- frobgev(data[,2], nsloc = nsloc2)
          mle.m1 <- mle.m1[c("loc","scale","shape")]
          mle.m2 <- mle.m2[c("loc","scale","shape")]
          qct[[i]] <- mtransform(qct[[i]], list(mle.m1, mle.m2), inv = TRUE)
        }
        else {
          if(!is.null(nsloc1) || !is.null(nsloc2))
            warning("ignoring `nsloc1' and `nsloc2' arguments")
          qct[[i]] <- mtransform(qct[[i]], kmar, inv = TRUE)
        }
      }
      qct[[i]][1,1] <- 1.5 * xlim[2]
      qct[[i]][nxv,2] <- 1.5 * ylim[2]
    }
    if((!is.null(nsloc1) || !is.null(nsloc2)) && !epmar && is.null(kmar)) {
        data <- fbvevd(data, model = "log", dep = 1, nsloc1 = nsloc1,
          nsloc2 = nsloc2, std.err = FALSE)$tdata
    }
    if(plot) {
      plot(data, xlab = xlab, ylab = ylab, xlim = xlim, ylim = ylim, ...)
      for(i in 1:np) lines(qct[[i]], lty = lty, lwd = lwd, col = col)
      return(invisible(qct))
    }
    if(add) {
      for(i in 1:np) lines(qct[[i]], lty = lty, lwd = lwd, col = col)
      return(invisible(qct))
    }
    qct
}

"amvnonpar"<- 
function(x = rep(1/3,3), data, epmar = FALSE, nsloc1 = NULL, nsloc2 = NULL,
    nsloc3 = NULL, method = c("pickands","deheuvels","halltajvidi"),
    madj = 0, kmar = NULL, plot = FALSE, col = heat.colors(12), blty = 0,
    grid = if(blty) 150 else 50, lower = 1/3, ord = 1:3, lab =
    as.character(1:3), lcex = 1)
{
    if(!plot) {
      if(is.vector(x)) x <- as.matrix(t(x))
      if(!is.matrix(x) || ncol(x) != 3)
        stop("`x' must be a vector/matrix with three elements/columns")
      if(any(x < 0, na.rm = TRUE))
        stop("`x' must be non-negative")
      rs <- rowSums(x)
      if(any(rs <= 0, na.rm = TRUE))
        stop("row(s) of `x' must have a positive sum")
      if(max(abs(rs[!is.na(rs)] - 1)) > 1e-6)
        warning("row(s) of `x' will be rescaled")
      x <- x/rs
    }
    if(missing(data) || ncol(data) != 3)
      stop("data must have three columns")
    if(epmar) {
      data <- apply(data, 2, rank, na.last = "keep")
      nasm <- apply(data, 2, function(x) sum(!is.na(x)))
      data <- data / rep(nasm+1, each = nrow(data))
      data <- -log(data)
    }
    else {
      if(is.null(kmar)) {
        if(!is.null(nsloc1)) {
            nsloc1 <- nsloc.transform(data, nsloc1)
            nslocmat1 <- cbind(1,as.matrix(nsloc1))
        }
        if(!is.null(nsloc2)) {
            nsloc2 <- nsloc.transform(data, nsloc2)
            nslocmat2 <- cbind(1,as.matrix(nsloc2))
        }
        if(!is.null(nsloc3)) {
            nsloc3 <- nsloc.transform(data, nsloc3)
            nslocmat3 <- cbind(1,as.matrix(nsloc3))
        }
        # Transform to exponential margins
        mle.m1 <- fgev(data[,1], nsloc = nsloc1, std.err = FALSE)$estimate
        loc.mle.m1 <- mle.m1[grep("^loc", names(mle.m1))]
        if(is.null(nsloc1)) loc.mle.m1 <- rep(loc.mle.m1, nrow(data))
        else loc.mle.m1 <- nslocmat1 %*% loc.mle.m1
        mle.m1 <- cbind(loc.mle.m1, mle.m1["scale"], mle.m1["shape"])
        
        mle.m2 <- fgev(data[,2], nsloc = nsloc2, std.err = FALSE)$estimate
        loc.mle.m2 <- mle.m2[grep("^loc", names(mle.m2))]
        if(is.null(nsloc2)) loc.mle.m2 <- rep(loc.mle.m2, nrow(data))
        else loc.mle.m2 <- nslocmat2 %*% loc.mle.m2
        mle.m2 <- cbind(loc.mle.m2, mle.m2["scale"], mle.m2["shape"])

        mle.m3 <- fgev(data[,3], nsloc = nsloc3, std.err = FALSE)$estimate
        loc.mle.m3 <- mle.m3[grep("^loc", names(mle.m3))]
        if(is.null(nsloc3)) loc.mle.m3 <- rep(loc.mle.m3, nrow(data))
        else loc.mle.m3 <- nslocmat3 %*% loc.mle.m3
        mle.m3 <- cbind(loc.mle.m3, mle.m3["scale"], mle.m3["shape"])

        data <- mtransform(data, list(mle.m1, mle.m2, mle.m3))
        # End transform
      }
      else {
        if(!is.null(nsloc1) || !is.null(nsloc2) || !is.null(nsloc3))
            warning("ignoring `nsloc1', `nsloc2' and `nsloc3' arguments")
        data <- mtransform(data, kmar)
      }
    }
    
    data <- na.omit(data)
    method <- match.arg(method)
    mpmin <- function(a,b,c) {
      a[a > b] <- b[a > b]
      a[a > c] <- c[a > c]
      a
    }
    
    if(method == "pickands" && (madj < 0.5)) {
      depfn <- function(x, data) {
        nn <- nrow(data)
        a <- numeric(nrow(x))
        for(i in 1:nrow(x))
          a[i] <- sum(mpmin(data[,1]/x[i,1], data[,2]/x[i,2], data[,3]/x[i,3]))
        a <- nn / a
        pmin(1, pmax(a, x[,1], x[,2], x[,3]))
      }
    }
    
    if((method == "deheuvels") ||
       ((method == "pickands") && (madj >= 0.5) && (madj < 1.5))) {
      depfn <- function(x, data) {
        nn <- nrow(data)
        a <- numeric(nrow(x))
        for(i in 1:nrow(x))
          a[i] <- sum(mpmin(data[,1]/x[i,1], data[,2]/x[i,2], data[,3]/x[i,3]))
        a <- nn / (a - x[,1] * sum(data[,1]) - x[,2] * sum(data[,2]) -
          x[,3] * sum(data[,3]) + nn)
        pmin(1, pmax(a, x[,1], x[,2], x[,3]))
      }
    }

    if(method == "halltajvidi" ||
       ((method == "pickands") && (madj >= 1.5))) {
      depfn <- function(x, data) {
        csum <- colSums(data)
        a <- numeric(nrow(x))
        for(i in 1:nrow(x))
          a[i] <- sum(mpmin(data[,1]/(csum[1] * x[i,1]),
            data[,2]/(csum[2] * x[i,2]), data[,3]/(csum[3] * x[i,3])))
        a <- 1 / a
        pmin(1, pmax(a, x[,1], x[,2], x[,3]))
      }
    }

    if(plot) {
      mz <- tvdepfn(depfn = depfn, col = col, blty = blty, grid = grid,
        lower = lower, ord = ord, lab = lab, lcex = lcex, data = data)
      return(invisible(mz))
    }
    depfn(x = x, data = data)
}

# Undocumented function; implements POT method
# See eqn (9.72) of Beirlant et al. (2004)
# Called from undocumented argument in abvnonpar
"abvnonparpot"<- 
function(x = 0.5, data, k = nrow(data)/4, epmar = FALSE, convex = FALSE,
         rev = FALSE, kmar = NULL, plot = FALSE, add = FALSE, lty = 1,
         lwd = 1, col = 1, blty = 3, blwd = 1, xlim = c(0,1), ylim =
         c(0.5,1), xlab = "t", ylab = "A(t)", ...)
{
    if(mode(x) != "numeric" || any(x < 0, na.rm=TRUE) ||
       any(x > 1, na.rm=TRUE)) stop("invalid argument for `x'")
    if(k >= nrow(data)) stop("k is too large")
    
    epdata <- apply(data, 2, rank, na.last = "keep")
    nasm <- apply(data, 2, function(x) sum(!is.na(x)))
    epdata <- epdata / rep(nasm+1, each = nrow(data))
    epdata <- -log(epdata)
    
    if(epmar) {data <- epdata}
    else {
      # Transform exceedances parametrically to exponential
      u1 <- sort(data[,1], decreasing = TRUE)[k+1]
      u2 <- sort(data[,2], decreasing = TRUE)[k+1]
      d1ab <- (data[,1] > u1); d2ab <- (data[,2] > u2)
      if(is.null(kmar)) {
        mle.m1 <- c(u1, fitted(fpot(data[d1ab,1], threshold = u1)))
        mle.m2 <- c(u2, fitted(fpot(data[d2ab,2], threshold = u2)))
        data[d1ab,1] <- mtransform(data[d1ab,1], mle.m1)
        data[d2ab,2] <- mtransform(data[d2ab,2], mle.m2)
      }
      else {
        data[d1ab,1] <- mtransform(data[d1ab,1], c(u1, kmar))
        data[d2ab,2] <- mtransform(data[d2ab,2], c(u2, kmar))
      }
      data[d1ab,1] <- -log(1 - k * data[d1ab,1] / nasm[1])
      data[d2ab,2] <- -log(1 - k * data[d2ab,2] / nasm[2])
      data[!d1ab, 1] <- epdata[!d1ab, 1]
      data[!d2ab, 2] <- epdata[!d2ab, 2]
      # End transform
    }
    
    if(rev) data <- data[,2:1]
    if(plot || add) x <- seq(0, 1, length = 100)
    mpmax <- function(a,b) {
      a[a < b] <- b[a < b]
      a
    }
    data <- na.omit(data)
    nn <- nrow(data)
    
    data <- 1/data
    rr <- rowSums(data)
    rrk <- sort(rr, decreasing = TRUE)[k+1]
    w1 <- data[,1]/rr

    if(!convex) {
      a <- numeric(length(x))
      for(i in 1:length(x))
        a[i] <- sum(mpmax(w1 * x[i], (1 - w1) * (1 - x[i]))[rr > rrk])
      a <- 2/k * a
      a <- pmin(1, pmax(a, x, 1-x))
    }
    else {
      # Not required: the method is always convex
      x2 <- seq(0, 1, length = 250)
      a <- numeric(250)
      for(i in 1:250)
        a[i] <- sum(mpmax(w1 * x2[i], (1 - w1) * (1 - x2[i]))[rr > rrk])
      a <- 2/k * a
      a <- pmin(1, pmax(a, x2, 1-x2))
      inch <- chull(x2, a)
      a <- a[inch] ; x2 <- x2[inch]
      a <- approx(x2, a, xout = x, method="linear")$y
    }
    if(plot || add) {
      bvdepfn(x, a, add, lty, lwd, col, blty, blwd, xlab, ylab, xlim, ylim, ...)  
      return(invisible(list(x = x, y = a)))
    } 
    a
}

