
"fextreme"<-
function(x, start, densfun, distnfun, ..., distn, mlen = 1, largest = TRUE,
         std.err = TRUE, corr = FALSE, method = "Nelder-Mead")
{
    if (missing(x) || length(x) == 0 || !is.numeric(x)) 
        stop("`x' must be a non-empty numeric object")
    if(any(is.na(x)))
        stop("`x' must not contain missing values")
    if (!is.list(start)) 
        stop("`start' must be a named list")
    call <- match.call()
    if(missing(densfun))
        densfun <- get(paste("d", distn, sep=""), mode="function")
    if(missing(distnfun))
        distnfun <- get(paste("p", distn, sep=""), mode="function")
    nllh <- function(p, ...) {
        dvec <- dens(p, ..., log = TRUE)
        if(any(is.infinite(dvec)))
            return(1e6)
        else 
            return(-sum(dvec))
    }
    nm <- names(start)
    l <- length(nm)
    f1 <- formals(densfun)
    f2 <- formals(distnfun)
    args <- names(f1)
    mtch <- match(nm, args)
    if (any(is.na(mtch))) 
        stop("`start' specifies unknown arguments")
    formals(densfun) <- c(f1[c(1, mtch)], f1[-c(1, mtch)])
    formals(distnfun) <- c(f2[c(1, mtch)], f2[-c(1, mtch)])
    dens <- function(p, x, densfun, distnfun, ...)
                dextreme(x, densfun, distnfun, p, ...)
    if(l > 1)
        body(dens) <- parse(text = paste("dextreme(x, densfun, distnfun,",
                            paste("p[",1:l,"]", collapse = ", "), ", ...)"))
    opt <- optim(start, nllh, x = x, hessian = TRUE, ...,
                 densfun = densfun, distnfun = distnfun, mlen = mlen,
                 largest = largest, method = method)
    if (opt$convergence != 0) {
        warning("optimization may not have succeeded")
        if(opt$convergence == 1) opt$convergence <- "iteration limit reached"
    }
    else opt$convergence <- "successful"
    if(std.err) {
        tol <- .Machine$double.eps^0.5
        var.cov <- qr(opt$hessian, tol = tol)
        if (var.cov$rank != ncol(var.cov$qr)) 
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
    structure(list(estimate = opt$par, std.err = std.err,
        deviance = 2*opt$value, corr = corr,
        convergence = opt$convergence, counts = opt$counts,
        message = opt$message, call = call, data = x,
        n = length(x)), class = c("extreme", "evd"))
}

"forder"<-
function(x, start, densfun, distnfun, ..., distn, mlen = 1, j = 1,
         largest = TRUE, std.err = TRUE, corr = FALSE, method = "Nelder-Mead")
{
    if (missing(x) || length(x) == 0 || !is.numeric(x)) 
        stop("`x' must be a non-empty numeric object")
    if(any(is.na(x)))
        stop("`x' must not contain missing values")
    if (!is.list(start)) 
        stop("`start' must be a named list")
    call <- match.call()
    if(missing(densfun))
        densfun <- get(paste("d", distn, sep=""), mode="function")
    if(missing(distnfun))
        distnfun <- get(paste("p", distn, sep=""), mode="function")
    nllh <- function(p, ...) {
        dvec <- dens(p, ..., log = TRUE)
        if(any(is.infinite(dvec)))
            return(1e6)
        else 
            return(-sum(dvec))
    }
    nm <- names(start)
    l <- length(nm)
    f1 <- formals(densfun)
    f2 <- formals(distnfun)
    args <- names(f1)
    mtch <- match(nm, args)
    if (any(is.na(mtch))) 
        stop("`start' specifies unknown arguments")
    formals(densfun) <- c(f1[c(1, mtch)], f1[-c(1, mtch)])
    formals(distnfun) <- c(f2[c(1, mtch)], f2[-c(1, mtch)])
    dens <- function(p, x, densfun, distnfun, ...)
                dorder(x, densfun, distnfun, p, ...)
    if(l > 1)
        body(dens) <- parse(text = paste("dorder(x, densfun, distnfun,",
                            paste("p[",1:l,"]", collapse = ", "), ", ...)"))
    opt <- optim(start, nllh, x = x, hessian = TRUE, ..., densfun = densfun,
                 distnfun = distnfun, mlen = mlen, j = j, largest = largest,
                 method = method)
    if (opt$convergence != 0) {
        warning("optimization may not have succeeded")
        if(opt$convergence == 1) opt$convergence <- "iteration limit reached"
    }
    else opt$convergence <- "successful"
    if(std.err) {
        tol <- .Machine$double.eps^0.5
        var.cov <- qr(opt$hessian, tol = tol)
        if (var.cov$rank != ncol(var.cov$qr)) 
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
    names(std.err) <- nm
    structure(list(estimate = opt$par, std.err = std.err,
        deviance = 2*opt$value, corr = corr,
        convergence = opt$convergence, counts = opt$counts,
        message = opt$message, call = call, data = x,
        n = length(x)), class = c("extreme", "evd"))
}

"fgev"<-
function(x, start, ..., nsloc = NULL, prob = NULL, std.err = TRUE,
    corr = FALSE, method = "BFGS", warn.inf = TRUE)
{
  call <- match.call()
  if(missing(x) || length(x) == 0 || !is.numeric(x)) 
    stop("`x' must be a non-empty numeric vector")
  if(is.null(prob)) {
    ft <- fgev.norm(x = x, start = start, ..., nsloc = nsloc, std.err =
      std.err, corr = corr, method = method, warn.inf = warn.inf)
  }
  else {
    if(length(prob) != 1 || !is.numeric(prob) || prob < 0 || prob > 1)
      stop("`prob' should be a probability in [0,1]")
    ft <- fgev.quantile(x = x, start = start, ..., nsloc = nsloc, prob = prob,
      std.err = std.err, corr = corr, method = method, warn.inf = warn.inf)
  }
  structure(c(ft, call = call), class = c("gev", "uvevd", "evd"))
}

"fgev.norm"<-
function(x, start, ..., nsloc = NULL, std.err = TRUE, corr = FALSE, method = "BFGS", warn.inf = TRUE)
{
    nlgev <- function(loc, scale, shape)
    { 
        if(scale <= 0) return(1e6)
        if(!is.null(nsloc)) {
            ns <- numeric(length(loc.param))
            for(i in 1:length(ns))
                ns[i] <- get(loc.param[i])
            loc <- drop(nslocmat %*% ns)
        }
        else loc <- rep(loc, length.out = length(x))
        .C("nlgev",
            x, n, loc, scale, shape, dns = double(1),
            PACKAGE = "evd")$dns
    }
    if(!is.null(nsloc)) {
        nsloc <- nsloc.transform(x, nsloc)
        nsloc <- nsloc[!is.na(x), ,drop = FALSE]
        nslocmat <- cbind(1,as.matrix(nsloc))
    }
    x <- as.double(x[!is.na(x)])
    n <- as.integer(length(x))
    loc.param <- paste("loc", c("",names(nsloc)), sep="")
    param <- c(loc.param, "scale", "shape")
    if(missing(start)) {
        start <- as.list(numeric(length(param)))
        names(start) <- param
        start$scale <- sqrt(6 * var(x))/pi
        start$loc <- mean(x) - 0.58 * start$scale
        start <- start[!(param %in% names(list(...)))]
    }
    if(!is.list(start)) 
        stop("`start' must be a named list")
    if(!length(start))
        stop("there are no parameters left to maximize over")
    nm <- names(start)
    l <- length(nm)
    f <- c(as.list(numeric(length(loc.param))), formals(nlgev)[2:3])
    names(f) <- param
    m <- match(nm, param)
    if(any(is.na(m))) 
        stop("`start' specifies unknown arguments")    
    formals(nlgev) <- c(f[m], f[-m])
    nllh <- function(p, ...) nlgev(p, ...)
    if(l > 1)
        body(nllh) <- parse(text = paste("nlgev(", paste("p[",1:l,
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
    if(!is.null(nsloc)) {
        trend <- param[paste("loc", names(nsloc), sep="")]
        trend <- drop(as.matrix(nsloc) %*% trend)
        x2 <- x - trend
    }
    else x2 <- x
    list(estimate = opt$par, std.err = std.err,
        fixed = unlist(fixed.param), param = param,
        deviance = 2*opt$value, corr = corr,
        convergence = opt$convergence, counts = opt$counts,
        message = opt$message,
        data = x, tdata = x2, nsloc = nsloc,
        n = length(x), prob = NULL, loc = param["loc"])
}

"fgev.quantile"<-
function(x, start, ..., nsloc = NULL, prob, std.err = TRUE, corr = FALSE, method = "BFGS", warn.inf = TRUE)
{
    nlgev <- function(quantile, scale, shape)
    {
        if(scale <= 0) return(1e6)
        quantile <- rep(quantile, length.out = length(x))
        if(prob == 0 && shape >= 0) return(1e6)
        if(prob == 1 && shape <= 0) return(1e6)
        if(shape == 0) loc <- quantile + scale * log(-log(1-prob))
        else loc <- quantile + scale/shape * (1 - (-log(1-prob))^(-shape))
        if(!is.null(nsloc)) {
            ns <- numeric(length(loc.param) - 1)
            for(i in 1:length(ns))
                ns[i] <- get(loc.param[i+1])
            loc <- drop(nslocmat %*% ns) + loc
        }
        if(any(is.infinite(loc))) return(1e6)
        .C("nlgev",
            x, n,
            loc, scale, shape, dns = double(1),
            PACKAGE = "evd")$dns
    }
    if(is.null(nsloc)) loc.param <- "quantile"
    else loc.param <- c("quantile", paste("loc", names(nsloc), sep=""))
    param <- c(loc.param, "scale", "shape")
    if(missing(start)) {
        start <- as.list(numeric(length(param)))
        names(start) <- param
        start$scale <- sqrt(6 * var(x, na.rm = TRUE))/pi
        start.loc <- mean(x, na.rm = TRUE) - 0.58 * start$scale
        start$quantile <- start.loc - start$scale * log(-log(1-prob))
        if(prob == 0) {
          fpft <- fgev(x = x, ..., nsloc = nsloc, prob = 0.001, std.err =
            std.err, corr = corr, method = method, warn.inf = warn.inf)
          start <- as.list(fitted(fpft))
        }
        if(prob == 1) {
          fpft <- fgev(x = x, ..., nsloc = nsloc, prob = 0.999, std.err =
            std.err, corr = corr, method = method, warn.inf = warn.inf)
          start <- as.list(fitted(fpft))
        }
        start <- start[!(param %in% names(list(...)))]
    }
    if(!is.list(start)) 
        stop("`start' must be a named list")
    if(!length(start))
        stop("there are no parameters left to maximize over")
    if(!is.null(nsloc)) {
        nsloc <- nsloc.transform(x, nsloc)
        nsloc <- nsloc[!is.na(x), ,drop = FALSE]
        nslocmat <- as.matrix(nsloc)
    }
    x <- as.double(x[!is.na(x)])
    n <- as.integer(length(x))
    nm <- names(start)
    l <- length(nm)
    f <- c(as.list(numeric(length(loc.param))), formals(nlgev)[2:3])
    names(f) <- param
    m <- match(nm, param)
    if(any(is.na(m))) 
        stop("`start' specifies unknown arguments")    
    formals(nlgev) <- c(f[m], f[-m])
    nllh <- function(p, ...) nlgev(p, ...)
    if(l > 1)
        body(nllh) <- parse(text = paste("nlgev(", paste("p[",1:l,
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
        .mat <- diag(1/std.err, nrow = length(std.err))
        if(corr) {
            corr <- structure(.mat %*% var.cov %*% .mat, dimnames = list(nm,nm))
            diag(corr) <- rep(1, length(std.err))
        }
        else corr <- NULL
    }
    else {
        std.err <- corr <- NULL
    }
    param <- c(opt$par, unlist(fixed.param))
    if(!is.null(nsloc)) {
        trend <- param[paste("loc", names(nsloc), sep="")]
        trend <- drop(as.matrix(nsloc) %*% trend)
        x2 <- x - trend
    }
    else x2 <- x
    if(param["shape"] == 0)
        loc <- param["quantile"] + param["scale"] * log(-log(1-prob))
    else
        loc <- param["quantile"] + param["scale"]/param["shape"] *
          (1 - (-log(1-prob))^(-param["shape"]))
    list(estimate = opt$par, std.err = std.err,
        fixed = unlist(fixed.param), param = param,
        deviance = 2*opt$value, corr = corr, convergence = opt$convergence,
        counts = opt$counts, message = opt$message,
        data = x, tdata = x2, nsloc = nsloc, n = length(x),
        prob = prob, loc = loc)
}

# Internal: robustly produces mles.
# Used for marginal parametric transforms.
"frobgev"<-
function(x, start, ..., nsloc = NULL) {
  mnx <- mean(x, na.rm = TRUE)
  sdx <- sqrt(var(x, na.rm = TRUE))
  mles <- fgev(x = (x-mnx)/sdx, start = start, ..., nsloc = nsloc,
    prob = NULL, std.err = FALSE, corr = FALSE)
  mles <- fitted(mles)
  mles["loc"] <- mnx + sdx * mles["loc"]
  mles["scale"] <- sdx * mles["scale"]
  trpar <- !(names(mles) %in% c("loc","scale","shape"))
  mles[trpar] <- sdx * mles[trpar]
  mles
}

"print.evd" <-  function(x, digits = max(3, getOption("digits") - 3), ...) 
{
    cat("\nCall:", deparse(x$call), "\n")
    cat("Deviance:", x$deviance, "\n")
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

"fitted.evd" <- function (object, ...) object$estimate
"std.errors" <- function (object, ...) UseMethod("std.errors")
"std.errors.evd" <- function (object, ...) object$std.err
"logLik.evd" <- function(object, ...) {
    val <- -deviance(object)/2
    attr(val, "df") <- length(fitted(object))
    class(val) <- "logLik"
    val
}

"nsloc.transform" <- 
function(x, nsloc)
{
    if(is.vector(nsloc))
        nsloc <- data.frame(trend = nsloc)
    if(!is.data.frame(nsloc))
        stop("`nsloc' must be a vector or data frame")
    if(is.null(dim(x))) ndat <- length(x)
    else ndat <- nrow(x)
    if(nrow(nsloc) != ndat)
        stop("`nsloc' and data are not compatible")
    nsloc
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





