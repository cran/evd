"pgpd" <-
function(q, loc = 0, scale = 1, shape = 0, lower.tail = TRUE)
{
    if(min(scale) <= 0) stop("invalid scale")
    if(length(shape) != 1) stop("invalid shape")
    q <- pmax(q - loc, 0)/scale
    if(shape == 0) p <- 1 - exp(-q)
    else {
        p <- pmax(1 + shape * q, 0)
        p <- 1 - p^(-1/shape)
    }
    if(!lower.tail) p <- 1 - p
    p
}

"dgpd"<-
function(x, loc = 0, scale = 1, shape = 0, log = FALSE)
{
    if(min(scale) <= 0) stop("invalid scale")
    if(length(shape) != 1) stop("invalid shape")
    d <- (x - loc)/scale
    nn <- length(d)
    scale <- rep(scale, length.out = nn)
    index <- (d > 0 & ((1 + shape * d) > 0)) | is.na(d)
    if(shape == 0) {
      d[index] <- log(1/scale[index]) - d[index]
      d[!index] <- -Inf
    }
    else {
        d[index] <- log(1/scale[index]) - (1/shape + 1) *
          log(1 + shape * d[index])
        d[!index] <- -Inf
    }
    if(!log) d <- exp(d)
    d
}

"qgpd"<-
function(p, loc = 0, scale = 1, shape = 0, lower.tail = TRUE)
{
    if(min(p, na.rm = TRUE) <= 0 || max(p, na.rm = TRUE) >=1)
        stop("`p' must contain probabilities in (0,1)")
    if(min(scale) < 0) stop("invalid scale")
    if(length(shape) != 1) stop("invalid shape")
    if(lower.tail) p <- 1 - p
    if(shape == 0) return(loc - scale*log(p))
    else return(loc + scale * (p^(-shape) - 1) / shape)
}

"rgpd"<-
function(n, loc = 0, scale = 1, shape = 0)
{
    if(min(scale) < 0) stop("invalid scale")
    if(length(shape) != 1) stop("invalid shape")
    if(shape == 0) return(loc + scale*rexp(n))
    else return(loc + scale * (runif(n)^(-shape) - 1) / shape)
}

"mrlplot"<-
function(data, tlim, nt = max(100, length(data)), lty = c(2,1,2), col = 1,
    conf = 0.95, main = "Mean Residual Life Plot", xlab = "Theshold",
    ylab = "Mean Excess", ...)
{
    data <- sort(data[!is.na(data)])
    nn <- length(data)
    if(nn <= 5) stop("`data' has too few non-missing values")
    if(missing(tlim)) {
      tlim <- c(data[1], data[nn - 4])
      tlim <- tlim - .Machine$double.eps^0.5
    }
    if(all(data <= tlim[2]))
      stop("upper limit for threshold is too high") 
    u <- seq(tlim[1], tlim[2], length = nt)
    x <- matrix(NA, nrow = nt, ncol = 3, dimnames = list(NULL,
        c("lower", "mrl", "upper")))
    for(i in 1:nt) {
        data <- data[data > u[i]]
	x[i,2] <- mean(data - u[i])
	sdev <- sqrt(var(data))
        sdev <- (qnorm((1 + conf)/2) * sdev)/sqrt(length(data))
	x[i,1] <- x[i,2] - sdev
	x[i,3] <- x[i,2] + sdev
    }
    matplot(u, x, type = "l", lty = lty, col = col, main = main,
            xlab = xlab, ylab = ylab, ...)
    invisible(list(x = u, y = x))
}

"tcplot"<-
function(data, tlim, model = c("gpd", "pp"), cmax = FALSE, r = 1, ulow =
    -Inf, rlow = 1, nt = 25, which = 1:npar, conf = 0.95, lty = 1, lwd = 1,
    type = "b", cilty = 1, ask = nb.fig < length(which) &&
    dev.interactive(), ...)
{
    model <- match.arg(model)
    u <- seq(tlim[1], tlim[2], length = nt)   
    locs <- scls <- shps <- matrix(NA, nrow = nt, ncol = 3)
    dimnames(locs) <- list(round(u,2), c("lower", "loc", "upper"))
    dimnames(shps) <- list(round(u,2), c("lower", "shape", "upper"))
    if(model == "gpd") {
      pname <- "mscale"
      npar <- 2
    }
    if(model == "pp") {
      pname <- "scale"
      npar <- 3
    }
    dimnames(scls) <- list(round(u,2), c("lower", pname, "upper"))
    z <- fpot(data, u[1], model = model, cmax = cmax, r = r, ulow = ulow,
              rlow = rlow, corr = TRUE, ...)
    stvals <- as.list(round(fitted(z), 3))
    for(i in 1:nt) {
      z <- fpot(data, u[i], model = model, start = stvals, cmax = cmax,
                r = r, ulow = ulow, rlow = rlow, corr = TRUE, ...)
      stvals <- as.list(fitted(z))
      mles <- fitted(z)
      stderrs <- std.errors(z)
      cnst <- qnorm((1 + conf)/2)
      shp <- mles["shape"]
      scl <- mles["scale"]
      shpse <- stderrs["shape"]
      sclse <- stderrs["scale"]
      if(model == "pp") {
        loc <- mles["loc"]  
        locse <- stderrs["loc"]
        locs[i,] <- c(loc - cnst*locse, loc, loc + cnst*locse)
      } 
      if(model == "gpd") {
        scl <- scl - shp*u[i]
        covar <- z$corr[1,2] * prod(stderrs)
        sclse <- sqrt(sclse^2 - 2*u[i]*covar + (u[i]*shpse)^2)
      }
      scls[i,] <- c(scl - cnst*sclse, scl, scl + cnst*sclse)
      shps[i,] <- c(shp - cnst*shpse, shp, shp + cnst*shpse)      
   }
   show <- rep(FALSE, npar)
   show[which] <- TRUE
   nb.fig <- prod(par("mfcol"))
   if (ask) {
     op <- par(ask = TRUE)
     on.exit(par(op))
   }
   if(model == "pp") {
     if(show[1]) {
       matplot(u, locs, type = "n", xlab = "Threshold",
               ylab = "Location")
       lines(u, locs[,2], lty = lty, lwd = lwd, type = type)
       segments(u, locs[,1], u, locs[,3], lty = cilty)
     }
     if(show[2]) {
       matplot(u, scls, type = "n", xlab = "Threshold",
               ylab = "Scale")
       lines(u, scls[,2], lty = lty, lwd = lwd, type = type)
       segments(u, scls[,1], u, scls[,3], lty = cilty)
     }
     if(show[3]) {
     matplot(u, shps, type = "n", xlab = "Threshold",
             ylab = "Shape")
     lines(u, shps[,2], lty = lty, lwd = lwd, type = type)
     segments(u, shps[,1], u, shps[,3], lty = cilty)
     }
     rtlist <- list(locs = locs, scales = scls, shapes = shps) 
   }
   if(model == "gpd") {
     if(show[1]) {
       matplot(u, scls, type = "n", xlab = "Threshold",
               ylab = "Modified Scale")
     lines(u, scls[,2], lty = lty, lwd = lwd, type = type)
     segments(u, scls[,1], u, scls[,3], lty = cilty)
     }
     if(show[2]) {
       matplot(u, shps, type = "n", xlab = "Threshold",
               ylab = "Shape")
       lines(u, shps[,2], lty = lty, lwd = lwd, type = type)
       segments(u, shps[,1], u, shps[,3], lty = cilty)
     }
     rtlist <- list(scales = scls, shapes = shps) 
   } 
   invisible(rtlist)
}

"fpot.norm"<-
function(x, threshold, model, start, npp = length(x), cmax = FALSE, r = 1, ulow = -Inf, rlow = 1, ..., std.err = TRUE, corr = FALSE, method = "BFGS", warn.inf = TRUE)
{
    if(model == "gpd") {
      nlpot <- function(scale, shape) { 
        .C("nlgpd",
            exceed, nhigh, threshold, scale, shape, dns = double(1),
            PACKAGE = "evd")$dns
      }
    }
    if(model == "pp") {
      nlpot <- function(loc, scale, shape) {
        .C("nlpp",
            exceed, nhigh, loc, scale, shape, threshold, nop,
            dns = double(1), PACKAGE = "evd")$dns
      }
    }
    nn <- length(x)
    nop <- as.double(nn/npp)
    if(cmax) {
      exceed <- clusters(x, u = threshold, r = r, ulow = ulow, rlow = rlow,
        cmax = TRUE, keep.names = FALSE)
      extind <- attributes(exceed)$acs
      exceed <- as.double(exceed)
      nhigh <- length(exceed) ; nat <- as.integer(nhigh * extind)
      extind <- 1/extind
    }
    else {
      extind <- r <- NULL
      high <- (x > threshold) & !is.na(x)
      exceed <- as.double(x[high])
      nhigh <- nat <- length(exceed)
    }
    if(!nhigh) stop("no data above threshold")
    pat <- nat/nn
    param <- c("scale", "shape")
    if(model == "pp") param <- c("loc", param)
    if(missing(start)) {
        if(model == "gpd") {
          start <- list(scale = 0, shape = 0)
          start$scale <- mean(exceed) - threshold
        }
        if(model == "pp") {
          start <- list(loc = 0, scale = 0, shape = 0)
          start$scale <- sqrt(6 * var(x))/pi
          start$loc <- mean(x) + (log(nop) - 0.58) * start$scale
        }
        start <- start[!(param %in% names(list(...)))]
    }
    if(!is.list(start)) 
        stop("`start' must be a named list")
    if(!length(start))
        stop("there are no parameters left to maximize over")
    nm <- names(start)
    l <- length(nm)
    f <- formals(nlpot)
    names(f) <- param
    m <- match(nm, param)
    if(any(is.na(m))) 
        stop("`start' specifies unknown arguments")    
    formals(nlpot) <- c(f[m], f[-m])
    nllh <- function(p, ...) nlpot(p, ...)
    if(l > 1)
        body(nllh) <- parse(text = paste("nlpot(", paste("p[",1:l,
            "]", collapse = ", "), ", ...)"))
    fixed.param <- list(...)[names(list(...)) %in% param]
    if(any(!(param %in% c(nm,names(fixed.param)))))
        stop("unspecified parameters")
    start.arg <- c(list(p = unlist(start)), fixed.param)
    if(warn.inf && do.call("nllh", start.arg) == 1e6)
        warning("negative log-likelihood is infinite at starting values")
    opt <- optim(start, nllh, hessian = TRUE, ..., method = method)
    if (opt$convergence != 0) {
        warning("optimization may not have succeeded")
        if(opt$convergence == 1) opt$convergence <- "iteration limit reached"
    }
    else opt$convergence <- "successful"
    if(std.err) {
        tol <- .Machine$double.eps^0.5
        var.cov <- qr(opt$hessian, tol = tol)
        if(var.cov$rank != ncol(var.cov$qr)) 
            stop("observed information matrix is singular; use std.err = FALSE")
        var.cov <- solve(var.cov, tol = tol)
        std.err <- diag(var.cov)
        if(any(std.err <= 0))
            stop("observed information matrix is singular; use std.err = FALSE")
        std.err <- sqrt(std.err)
        names(std.err) <- nm
        if(corr) {
            .mat <- diag(1/std.err, nrow = length(std.err))
            corr <- structure(.mat %*% var.cov %*% .mat, dimnames = list(nm,nm))
            diag(corr) <- rep(1, length(std.err))
        }
        else corr <- NULL
    }
    else std.err <- corr <- NULL
    param <- c(opt$par, unlist(fixed.param))
   if(model == "gpd") scale <- param["scale"]
   if(model == "pp") scale <- param["scale"] + param["shape"] * (threshold -
     param["loc"])
    
    list(estimate = opt$par, std.err = std.err, fixed = unlist(fixed.param),
        param = param, deviance = 2*opt$value, corr = corr, convergence =
        opt$convergence, counts = opt$counts, message = opt$message,
        threshold = threshold, r = r, ulow = ulow, rlow = rlow, npp = npp,
        nhigh = nhigh, nat = nat, pat = pat, extind = extind,
        data = x, exceedances = exceed, mper = NULL, scale = scale)
}

"fpot.quantile"<-
function(x, threshold, start, npp = length(x), cmax = FALSE, r = 1, ulow = -Inf, rlow = 1, mper, ..., std.err = TRUE, corr = FALSE, method = "BFGS", warn.inf = TRUE)
{
    nlpot <- function(rlevel, shape)
    {
        if(is.infinite(mper) && shape >= 0) return(1e6)
        rlevel <- rlevel - threshold
        if(shape == 0) scale <- rlevel / log(adjmper)
        else scale <- shape * rlevel / (adjmper^shape - 1)
        .C("nlgpd",
            exceed, nhigh, threshold, scale, shape, dns = double(1),
            PACKAGE = "evd")$dns
    }
    nn <- length(x)
    if(cmax) {
      exceed <- clusters(x, u = threshold, r = r, ulow = ulow, rlow = rlow,
        cmax = TRUE, keep.names = FALSE)
      extind <- attributes(exceed)$acs
      exceed <- as.double(exceed)
      nhigh <- length(exceed) ; nat <- as.integer(nhigh * extind)
      extind <- 1/extind
    }
    else {
      extind <- r <- NULL
      high <- (x > threshold) & !is.na(x)
      exceed <- as.double(x[high])
      nhigh <- nat <- length(exceed)
    }
    if(!nhigh) stop("no data above threshold")
    pat <- nat/nn
    adjmper <- mper * npp * nhigh/nn
    if(adjmper <= 1) stop("`mper' is too small")
    param <- c("rlevel", "shape")
    if(missing(start)) {
        start <- list(rlevel = 0, shape = 0)
        stscale <- mean(exceed) - threshold
        start$rlevel <- threshold + stscale*log(adjmper)
        if(is.infinite(mper)) {
          stmp <- 100/(npp * nhigh/nn)
          fpft <- fpot(x = x, threshold = threshold, npp = npp, cmax =
            cmax, r = r, ulow = ulow, rlow = rlow, mper = stmp, ...,
            std.err = std.err, corr = corr, method = method, warn.inf =
            warn.inf)
          start <- as.list(fitted(fpft))
        }
        start <- start[!(param %in% names(list(...)))]
    }
    if(!is.list(start)) 
        stop("`start' must be a named list")
    if(!length(start))
        stop("there are no parameters left to maximize over")
    nm <- names(start)
    l <- length(nm)
    f <- formals(nlpot)
    names(f) <- param
    m <- match(nm, param)
    if(any(is.na(m))) 
        stop("`start' specifies unknown arguments")    
    formals(nlpot) <- c(f[m], f[-m])
    nllh <- function(p, ...) nlpot(p, ...)
    if(l > 1)
        body(nllh) <- parse(text = paste("nlpot(", paste("p[",1:l,
            "]", collapse = ", "), ", ...)"))
    fixed.param <- list(...)[names(list(...)) %in% param]
    if(any(!(param %in% c(nm,names(fixed.param)))))
        stop("unspecified parameters")
    start.arg <- c(list(p = unlist(start)), fixed.param)
    if(warn.inf && do.call("nllh", start.arg) == 1e6)
        warning("negative log-likelihood is infinite at starting values")
    opt <- optim(start, nllh, hessian = TRUE, ..., method = method)
    if (opt$convergence != 0) {
        warning("optimization may not have succeeded")
        if(opt$convergence == 1) opt$convergence <- "iteration limit reached"
    }
    else opt$convergence <- "successful"
    if(std.err) {
        tol <- .Machine$double.eps^0.5
        var.cov <- qr(opt$hessian, tol = tol)
        if(var.cov$rank != ncol(var.cov$qr)) 
            stop("observed information matrix is singular; use std.err = FALSE")
        var.cov <- solve(var.cov, tol = tol)
        std.err <- diag(var.cov)
        if(any(std.err <= 0))
            stop("observed information matrix is singular; use std.err = FALSE")
        std.err <- sqrt(std.err)
        names(std.err) <- nm
        if(corr) {
            .mat <- diag(1/std.err, nrow = length(std.err))
            corr <- structure(.mat %*% var.cov %*% .mat, dimnames = list(nm,nm))
            diag(corr) <- rep(1, length(std.err))
        }
        else corr <- NULL
    }
    else std.err <- corr <- NULL
    param <- c(opt$par, unlist(fixed.param))
    rlevel <- param["rlevel"] - threshold
    if(param["shape"] == 0) scale <- rlevel / log(adjmper)
    else scale <- param["shape"] * rlevel / (adjmper^param["shape"] - 1) 
    list(estimate = opt$par, std.err = std.err, fixed = unlist(fixed.param),
        param = param, deviance = 2*opt$value, corr = corr, convergence =
        opt$convergence, counts = opt$counts, message = opt$message,
        threshold = threshold, r = r, ulow = ulow, rlow = rlow, npp = npp,
        nhigh = nhigh, nat = nat, pat = pat, extind = extind,
        data = x, exceedances = exceed, mper = mper, scale = scale)
}

"fpot"<-
function(x, threshold, model = c("gpd", "pp"), start, npp = length(x), cmax = FALSE, r = 1, ulow = -Inf, rlow = 1, mper = NULL, ..., std.err = TRUE, corr = FALSE, method = "BFGS", warn.inf = TRUE)
{
  call <- match.call()
  model <- match.arg(model)
  if(missing(x) || length(x) == 0 || mode(x) != "numeric") 
    stop("`x' must be a non-empty numeric vector")
  if(missing(threshold) || length(threshold) != 1 ||
     mode(threshold) != "numeric") 
    stop("`threshold' must be a numeric value")
  threshold <- as.double(threshold)
  if(is.null(mper)) {
    ft <- fpot.norm(x = x, threshold = threshold, model = model, start = start,
      npp = npp, cmax = cmax, r = r, ulow = ulow, rlow = rlow, ...,
      std.err = std.err, corr = corr, method = method, warn.inf = warn.inf)
  }
  else {
    if(model == "pp")
      stop("`mper' cannot be specified in point process models")
    ft <- fpot.quantile(x = x, threshold = threshold, start =
      start, npp = npp, cmax = cmax, r = r, ulow = ulow, rlow = rlow, ...,
      mper = mper, std.err = std.err, corr = corr, method = method,
      warn.inf = warn.inf)
  }
  structure(c(ft, call = call), class = c("pot", "uvevd", "evd"))
}

"print.pot" <-  function(x, digits = max(3, getOption("digits") - 3), ...) 
{
    cat("\nCall:", deparse(x$call), "\n")
    cat("Deviance:", x$deviance, "\n")

    cat("\nThreshold:", round(x$threshold, digits), "\n")
    cat("Number Above:", x$nat, "\n")
    cat("Proportion Above:", round(x$pat, digits), "\n")
    if(!is.null(x$extind)) {
      cat("\nClustering Interval:", x$r, "\n")
      if(is.finite(x$ulow)) {
        cat("Lower Threshold:", round(x$ulow, digits), "\n")
        cat("Lower Clustering Interval:", x$rlow, "\n")
      }
      cat("Number of Clusters:", x$nhigh, "\n")
      cat("Extremal Index:", round(x$extind, digits), "\n")
    }
    
    cat("\nEstimates\n") 
    print.default(format(x$estimate, digits = digits), print.gap = 2, 
        quote = FALSE)
    if(!is.null(x$std.err)) {
    cat("\nStandard Errors\n")
    print.default(format(x$std.err, digits = digits), print.gap = 2, 
        quote = FALSE)
    }
    if(!is.null(x$corr)) {
    cat("\nCorrelations\n")
    print.default(format(x$corr, digits = digits), print.gap = 2, 
        quote = FALSE)
    }
    cat("\nOptimization Information\n")
    cat("  Convergence:", x$convergence, "\n")
    cat("  Function Evaluations:", x$counts["function"], "\n")
    if(!is.na(x$counts["gradient"]))
        cat("  Gradient Evaluations:", x$counts["gradient"], "\n")
    if(!is.null(x$message)) cat("  Message:", x$message, "\n")
    cat("\n")
    invisible(x)
}

"qq.pot" <-  function(x, ci = TRUE, main = "Quantile Plot", xlab = "Model", ylab = "Empirical", ...)
{
    quant <- qgpd(ppoints(x$nhigh), loc = x$threshold,
                 scale = x$scale, shape = x$param["shape"])
    if(!ci) {
      plot(quant, sort(x$exceedances), main = main, xlab = xlab,
           ylab = ylab, ...)
      abline(0, 1)
    }
    else {
      samp <- rgpd(x$nhigh*99, loc = x$threshold,
                 scale = x$scale, shape = x$param["shape"])
      samp <- matrix(samp, x$nhigh, 99)
      samp <- apply(samp, 2, sort)
      samp <- apply(samp, 1, sort)
      env <- t(samp[c(3,97),])
      rs <- sort(x$exceedances)
      matplot(quant, cbind(rs,env), main = main, xlab = xlab, ylab = ylab,
              type = "pnn", pch = 4, ...)
      xyuser <- par("usr")
      smidge <- min(diff(c(xyuser[1], quant, xyuser[2])))/2
      segments(quant-smidge, env[,1], quant+smidge, env[,1])
      segments(quant-smidge, env[,2], quant+smidge, env[,2])
      abline(0, 1)
    }
    invisible(list(x = quant, y = sort(x$exceedances)))
}

"pp.pot" <-  function(x, ci = TRUE, main = "Probability Plot", xlab = "Empirical", ylab = "Model", ...)
{
    ppx <- ppoints(x$nhigh)
    probs <- pgpd(sort(x$exceedances), loc = x$threshold,
                 scale = x$scale, shape = x$param["shape"])
    if(!ci) {
        plot(ppx, probs, main = main, xlab = xlab, ylab = ylab, ...)
        abline(0, 1)
    }
    else {
        samp <- rgpd(x$nhigh*99, loc = x$threshold,
                   scale = x$scale, shape = x$param["shape"])
        samp <- matrix(samp, x$nhigh, 99)
        samp <- apply(samp, 2, sort)
        samp <- apply(samp, 1, sort)
        env <- t(samp[c(3,97),])
        env[,1] <- pgpd(env[,1], loc = x$threshold,
                    scale = x$scale, shape = x$param["shape"])
        env[,2] <- pgpd(env[,2], loc = x$threshold,
                    scale = x$scale, shape = x$param["shape"])
        matplot(ppx, cbind(probs, env), main = main, xlab = xlab,
                ylab = ylab, type = "pnn", pch = 4, ...)
        xyuser <- par("usr")
        smidge <- min(diff(c(xyuser[1], ppx, xyuser[2])))/2
        segments(ppx-smidge, env[,1], ppx+smidge, env[,1])
        segments(ppx-smidge, env[,2], ppx+smidge, env[,2])
        abline(0, 1)
    }
    invisible(list(x = ppoints(x$nhigh), y = probs))
}

"rl.pot" <- function(x, ci = TRUE, main = "Return Level Plot", xlab = "Return Period", ylab = "Return Level", ...)
{
    rpstmfc <- c(1.001,10^(seq(0,3,len=200))[-1])
    rlev <- qgpd(1/rpstmfc, loc = x$threshold, scale = x$scale,
              shape = x$param["shape"], lower.tail = FALSE)
    mfc <- x$npp * x$nhigh/length(x$data)
    rps <- rpstmfc/mfc
    ppx <- 1/(mfc * (1 - ppoints(x$nhigh)))
    if(!ci) {
        plot(ppx, sort(x$exceedances), log = "x", main =
             main, xlab = xlab, ylab = ylab, ...)
        lines(rps, rlev)
    }
    else {
        samp <- rgpd(x$nhigh*99, loc = x$threshold,
                   scale = x$scale, shape = x$param["shape"])
        samp <- matrix(samp, x$nhigh, 99)
        samp <- apply(samp, 2, sort)
        samp <- apply(samp, 1, sort)
        env <- t(samp[c(3,97),])
        rs <- sort(x$exceedances)
        matplot(ppx, cbind(rs,env), main = main, xlab = xlab,
                ylab = ylab, type = "pnn", pch = 4, log = "x", ...)
        lines(rps, rlev)
        xyuser <- par("usr")
        smidge <- min(diff(c(xyuser[1], log10(ppx), xyuser[2])))/2
        segments(ppx*exp(-smidge), env[,1], ppx*exp(smidge), env[,1])
        segments(ppx*exp(-smidge), env[,2], ppx*exp(smidge), env[,2])
    }
    invisible(list(x = rps, y = rlev))
}

"dens.pot" <-  function(x, adjust = 1, nplty = 2, jitter = FALSE, main = "Density Plot", xlab = "Quantile", ylab = "Density", ...)
{
    xlimit <- range(x$exceedances)
    xlimit[2] <- xlimit[2] + diff(xlimit) * 0.075
    xvec <- seq(xlimit[1], xlimit[2], length = 100)
    dens <- dgpd(xvec, loc = x$threshold, scale = x$scale,
                shape = x$param["shape"])
    plot(spline(xvec, dens), main = main, xlab = xlab, ylab = ylab,
         type = "l", ...)
    if(jitter) rug(jitter(x$exceedances))
    else rug(x$exceedances)
    flipexceed <- c(x$exceedances, 2*x$threshold - x$exceedances)
    flipdens <- density(flipexceed, adjust = adjust, from = xlimit[1],
        to = xlimit[2])
    flipdens$y <- 2*flipdens$y
    lines(flipdens, lty = nplty)
    invisible(list(x = xvec, y = dens))
}

"clusters"<-
function(data, u, r = 1, ulow = -Inf, rlow = 1, cmax = FALSE, keep.names = TRUE,
    plot = FALSE, xdata = seq(along = data), lvals = TRUE, lty = 1, lwd = 1,
    pch = par("pch"), col = if(n > 250) NULL else "grey", xlab = "Index",
    ylab = "Data", ...)
{
    n <- length(data)
    if(length(u) != 1) u <- rep(u, length.out = n)
    if(length(ulow) != 1) ulow <- rep(ulow, length.out = n)
    if(any(ulow > u)) stop("`u' cannot be less than `ulow'")
    if(is.null(names(data)) && keep.names) names(data) <- 1:n
    if(!keep.names) names(data) <- NULL
    high <- as.double((data > u) & !is.na(data))
    high2 <- as.double((data > ulow) | is.na(data))
    clstrs <- .C("clusters", high, high2, n, as.integer(r),
        as.integer(rlow), clstrs = double(3*n), PACKAGE = "evd")$clstrs
    clstrs <- matrix(clstrs, n, 3)
    start <- clstrs[,2] ; end <- clstrs[,3]
    splvec <- clstrs[,1]
    start <- as.logical(start)
    end <- as.logical(end)
    clstrs <- split(data, splvec)
    names(clstrs) <- paste("cluster", names(clstrs), sep = "")
    if(any(!splvec)) clstrs <- clstrs[-1]
    nclust <- length(clstrs)
    acs <- sum(high)/nclust
    if(plot) {
      if(length(xdata) != length(data))
        stop("`xdata' and `data' have different lengths")
      if(any(is.na(xdata)))
        stop("`xdata' cannot contain missing values")
      if(any(duplicated(xdata)))
        stop("`xdata' cannot contain duplicated values")
      eps <- min(diff(xdata))/2
      start <- xdata[start] - eps
      end <- xdata[end] + eps
      plot(xdata, data, xlab = xlab, ylab = ylab, type = "n", ...)
      if(!is.null(col)) {
        for(i in 1:nclust) {
          xvl <- c(start[i], end[i], end[i], start[i])
          polygon(xvl, rep(par("usr")[3:4], each = 2), col = col)
        }
      }
      if(length(u) == 1) abline(h = u, lty = lty, lwd = lwd)
      else lines(xdata, u, lty = lty, lwd = lwd)
      if(lvals) {
        if(length(ulow) == 1) abline(h = ulow, lty = lty, lwd = lwd)
        else lines(xdata, ulow, lty = lty, lwd = lwd)
      }
      else {
        high <- as.logical(high)
        xdata <- xdata[high]
        data <- data[high]
      }
      points(xdata, data, pch = pch)
    }
    if(cmax) {
      if(keep.names)
        nmcl <- unlist(lapply(clstrs, function(x) names(x)[which.max(x)]))
      clstrs <- as.numeric(unlist(lapply(clstrs, max, na.rm = TRUE)))
      if(keep.names) names(clstrs) <- nmcl 
    }
    attributes(clstrs)$acs <- acs
    if(plot) return(invisible(clstrs))   
    clstrs
}

"exi"<-
function(data, u, r = 1, ulow = rep(-Inf, ncol(u)), rlow = rep(1, length(r)),
         dimnames = list(NULL, NULL), drop = TRUE)
{
    if(!is.matrix(u)) u <- t(as.matrix(u))
    if(!is.matrix(ulow)) ulow <- t(as.matrix(ulow))
    if(ncol(u) != ncol(ulow)) stop("`u' and `ulow' are not compatible")
    if(length(r) != length(rlow)) stop("`r' and `rlow' are not compatible")
    n <- ncol(u) ; m <- length(r)
    mat <- matrix(0, n, m, dimnames = dimnames)
    for(i in 1:n) {
       for(j in 1:m) {
          clstrs <- clusters(data, u = u[,i], r = r[j], ulow = ulow[,i],
            rlow = rlow[j], keep.names = FALSE)
          mat[i, j] <- attributes(clstrs)$acs
       }
    }
    if(drop) mat <- drop(mat)
    1/mat
}

