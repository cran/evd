"rfrechet"<-
function(n, loc = 0, scale = 1, shape = 1)
{
    if(min(scale) < 0 || min(shape) <= 0) stop("invalid arguments")
    loc + scale * rexp(n)^(-1/shape)
}

"rgumbel"<-
function(n, loc = 0, scale = 1)
{
    rgev(n, loc = loc, scale = scale, shape = 0)
}

"rrweibull"<-
function(n, loc = 0, scale = 1, shape = 1)
{
    if(min(scale) < 0 || min(shape) <= 0) stop("invalid arguments")
    loc - scale * rexp(n)^(1/shape)
}

"rgev"<-
function(n, loc = 0, scale = 1, shape = 0)
{
    if(min(scale) < 0) stop("invalid scale")
    if(length(shape) != 1) stop("invalid shape")
    if(shape == 0) return(loc - scale * log(rexp(n)))
    else return(loc + scale * (rexp(n)^(-shape) - 1)/shape)
}

"rorder"<-
function(n, quantfun, ..., distn,  mlen = 1, j = 1, largest = TRUE)
{
    if(mode(mlen) != "numeric" || length(mlen) != 1 || mlen < 1 ||
       mlen %% 1 != 0) 
        stop("`mlen' must be a non-negative integer")
    if(mode(j) != "numeric" || length(j) != 1 || j < 1 || j %% 1 != 0) 
        stop("`j' must be a non-negative integer")
    if(j > mlen)
        stop("`j' cannot be greater than `mlen'")
    if(!largest) j <- mlen+1-j
    if(missing(quantfun))
        quantfun <- get(paste("q", distn, sep=""), mode="function")
    quantfun(rbeta(n, mlen+1-j, j), ...)
}

"rext"<-
function(n, quantfun, ..., distn, mlen = 1, largest = TRUE)
{
    if(mode(mlen) != "numeric" || length(mlen) != 1 || mlen < 1 ||
       mlen %% 1 != 0) 
        stop("`mlen' must be a non-negative integer")
    if(missing(quantfun))
        quantfun <- get(paste("q", distn, sep=""), mode="function")
    if(largest)
        quantfun(rbeta(n, mlen, 1), ...)
    else
        quantfun(rbeta(n, 1, mlen), ...) 
}

"qfrechet"<-
function(p, loc = 0, scale = 1, shape = 1, lower.tail = TRUE)
{
    if(min(p, na.rm = TRUE) <= 0 || max(p, na.rm = TRUE) >=1)
        stop("`p' must contain probabilities in (0,1)")
    if(min(scale) < 0 || min(shape) <= 0) stop("invalid arguments")
    if(!lower.tail) p <- 1 - p
    loc + scale * (-log(p))^(-1/shape)
}

"qgumbel"<-
function(p, loc = 0, scale = 1, lower.tail = TRUE)
{
    qgev(p, loc = loc, scale = scale, shape = 0, lower.tail = lower.tail)
}

"qrweibull"<-
function(p, loc = 0, scale = 1, shape = 1, lower.tail = TRUE)
{
    if(min(p, na.rm = TRUE) <= 0 || max(p, na.rm = TRUE) >=1)
        stop("`p' must contain probabilities in (0,1)")
    if(min(scale) < 0 || min(shape) <= 0) stop("invalid arguments")
    if(!lower.tail) p <- 1 - p
    loc - scale * (-log(p))^(1/shape)
}

"qgev"<-
function(p, loc = 0, scale = 1, shape = 0, lower.tail = TRUE)
{
    if(min(p, na.rm = TRUE) <= 0 || max(p, na.rm = TRUE) >=1)
        stop("`p' must contain probabilities in (0,1)")
    if(min(scale) < 0) stop("invalid scale")
    if(length(shape) != 1) stop("invalid shape")
    if(!lower.tail) p <- 1 - p
    if(shape == 0) return(loc - scale * log(-log(p)))
    else return(loc + scale * ((-log(p))^(-shape) - 1)/shape)
}

"qext"<-
function(p, quantfun, ..., distn, mlen = 1, largest = TRUE, lower.tail = TRUE)
{
    if(min(p, na.rm = TRUE) <= 0 || max(p, na.rm = TRUE) >=1)
        stop("`p' must contain probabilities in (0,1)")
    if(mode(mlen) != "numeric" || length(mlen) != 1 || mlen < 1 ||
       mlen %% 1 != 0) 
        stop("`mlen' must be a non-negative integer")
    if(missing(quantfun))
        quantfun <- get(paste("q", distn, sep=""), mode="function")
    if(!lower.tail) p <- 1 - p
    if(largest) 
        quantfun(p^(1/mlen), ...)
    else
        quantfun(1-(1-p)^(1/mlen), ...)
}

"pfrechet"<-
function(q, loc = 0, scale = 1, shape = 1, lower.tail = TRUE)
{
    if(min(scale) <= 0 || min(shape) <= 0) stop("invalid arguments")
    q <- pmax((q - loc)/scale,0)
    p <- exp(-q^(-shape))
    if(!lower.tail) p <- 1 - p
    p
}

"pgumbel"<-
function(q, loc = 0, scale = 1, lower.tail = TRUE)
{
    pgev(q, loc = loc, scale = scale, shape = 0, lower.tail = lower.tail)
}

"prweibull"<-
function(q, loc = 0, scale = 1, shape = 1, lower.tail = TRUE)
{
    if(min(scale) <= 0 || min(shape) <= 0) stop("invalid arguments")
    q <- pmin((q - loc)/scale,0)
    p <- exp(-(-q)^shape)
    if(!lower.tail) p <- 1 - p
    p
}

"pgev"<-
function(q, loc = 0, scale = 1, shape = 0, lower.tail = TRUE)
{
    if(min(scale) <= 0) stop("invalid scale")
    if(length(shape) != 1) stop("invalid shape")
    q <- (q - loc)/scale
    if(shape == 0) p <- exp(-exp(-q))
    else p <- exp( - pmax((1 + shape * q),0)^(-1/shape))
    if(!lower.tail) p <- 1 - p
    p
}

"porder"<-
function(q, distnfun, ..., distn, mlen = 1, j = 1, largest = TRUE,
         lower.tail = TRUE)
{
    if(mode(mlen) != "numeric" || length(mlen) != 1 || mlen < 1 ||
       mlen %% 1 != 0) 
        stop("`mlen' must be a non-negative integer")
    if(mode(j) != "numeric" || length(j) != 1 || j < 1 || j %% 1 != 0) 
        stop("`j' must be a non-negative integer")
    if(j > mlen)
        stop("`j' cannot be greater than `mlen'")
    lachooseb <- function(a,b) lgamma(a+1) - lgamma(b+1) - lgamma(a-b+1)
    if(largest) svec <- (mlen+1-j):mlen
    else  svec <- 0:(j-1)
    if(missing(distnfun))
        distnfun <- get(paste("p", distn, sep=""), mode="function")
    distn <- distnfun(q, ...)
    store <- matrix(0,nrow=length(q),ncol=j)
    for(k in 1:j)
        store[,k] <- exp(lachooseb(mlen,svec[k]) + svec[k]*log(distn) +
                         (mlen-svec[k])*log(1-distn))
    p <- apply(store,1,sum)
    if(largest != lower.tail) p <- 1 - p
    p
}

"pext"<-
function(q, distnfun, ..., distn, mlen = 1, largest = TRUE, lower.tail = TRUE)
{
    if(mode(mlen) != "numeric" || length(mlen) != 1 || mlen < 1 ||
       mlen %% 1 != 0) 
        stop("`mlen' must be a non-negative integer")
    if(missing(distnfun))
        distnfun <- get(paste("p", distn, sep=""), mode="function")
    distn <- distnfun(q, ...)
    if(!largest) distn <- 1-distn
    p <- distn^mlen
    if(largest != lower.tail) p <- 1 - p
    p
}

"dfrechet"<-
function(x, loc = 0, scale = 1, shape = 1, log = FALSE)
{
    if(min(scale) <= 0 || min(shape) <= 0) stop("invalid arguments")
    x <- (x - loc)/scale
    xpos <- x[x>0 | is.na(x)]
    nn <- length(x)
    scale <- rep(scale, length.out = nn)[x>0 | is.na(x)]
    shape <- rep(shape, length.out = nn)[x>0 | is.na(x)]
    d <- numeric(nn)
    d[x>0 | is.na(x)] <- log(shape/scale) - (1+shape) * log(xpos) -
         xpos^(-shape)
    d[x<=0 & !is.na(x)] <- -Inf
    if(!log) d <- exp(d)
    d
}

"dgumbel"<-
function(x, loc = 0, scale = 1, log = FALSE)
{
    dgev(x, loc = loc, scale = scale, shape = 0, log = log)
}

"drweibull"<-
function(x, loc = 0, scale = 1, shape = 1, log = FALSE)
{
    if(min(scale) <= 0 || min(shape) <= 0) stop("invalid arguments")
    x <- (x - loc)/scale
    xneg <- x[x<0 | is.na(x)]
    nn <- length(x)
    scale <- rep(scale, length.out = nn)[x<0 | is.na(x)]
    shape <- rep(shape, length.out = nn)[x<0 | is.na(x)]
    d <- numeric(nn)
    d[x<0 | is.na(x)] <- log(shape/scale) + (shape-1) * log(-xneg) -
        (-xneg)^shape
    d[x>=0 & !is.na(x)] <- -Inf
    if(!log) d <- exp(d)
    d
}

"dgev"<-
function(x, loc = 0, scale = 1, shape = 0, log = FALSE)
{
    if(min(scale) <= 0) stop("invalid scale")
    if(length(shape) != 1) stop("invalid shape")
    x <- (x - loc)/scale
    if(shape == 0)
        d <- log(1/scale) - x - exp(-x) 
    else {
        nn <- length(x)
        xx <- 1 + shape*x
        xxpos <- xx[xx>0 | is.na(xx)]
        scale <- rep(scale, length.out = nn)[xx>0 | is.na(xx)]
        d <- numeric(nn)
        d[xx>0 | is.na(xx)] <- log(1/scale) - xxpos^(-1/shape) -
            (1/shape + 1)*log(xxpos)
        d[xx<=0 & !is.na(xx)] <- -Inf
    }  
    if(!log) d <- exp(d)
    d
}

"dorder"<-
function(x, densfun, distnfun, ..., distn, mlen = 1, j = 1, largest = TRUE,
         log = FALSE)
{
    if(mode(mlen) != "numeric" || length(mlen) != 1 || mlen < 1 ||
       mlen %% 1 != 0) 
        stop("`mlen' must be a non-negative integer")
    if(mode(j) != "numeric" || length(j) != 1 || j < 1 || j %% 1 != 0) 
        stop("`j' must be a non-negative integer")
    if(j > mlen)
        stop("`j' cannot be greater than `mlen'")
    if(!largest) j <- mlen + 1 - j
    if(missing(densfun))
        densfun <- get(paste("d", distn, sep=""), mode="function")
    if(missing(distnfun))
        distnfun <- get(paste("p", distn, sep=""), mode="function")
    dens <- densfun(x, ..., log = TRUE)
    distn <- distnfun(x, ...)[!is.infinite(dens)]
    distn <- (mlen-j) * log(distn) + (j-1) * log(1-distn)
    comb <- lgamma(mlen+1) - lgamma(j) - lgamma(mlen-j+1)
    d <- numeric(length(x))
    d[!is.infinite(dens)] <- comb + dens[!is.infinite(dens)] + distn
    d[is.infinite(dens)] <- -Inf
    if(!log) d <- exp(d)
    d
}

"dext"<-
function(x, densfun, distnfun, ..., distn, mlen = 1, largest = TRUE, log = FALSE)
{
    if(mode(mlen) != "numeric" || length(mlen) != 1 || mlen < 1 ||
       mlen %% 1 != 0) 
        stop("`mlen' must be a non-negative integer")
    if(!largest) j <- mlen + 1 - j
    if(missing(densfun))
        densfun <- get(paste("d", distn, sep=""), mode="function")
    if(missing(distnfun))
        distnfun <- get(paste("p", distn, sep=""), mode="function")
    dens <- densfun(x, ..., log = TRUE)
    distn <- distnfun(x, ...)[!is.infinite(dens)]
    if(!largest) distn <- 1 - distn
    distn <- (mlen-1) * log(distn)
    d <- numeric(length(x))
    d[!is.infinite(dens)] <- log(mlen) + dens[!is.infinite(dens)] + distn
    d[is.infinite(dens)] <- -Inf
    if(!log) d <- exp(d)
    d
}

"fext"<-
function(x, start, densfun, distnfun, ..., distn, mlen = 1, largest = TRUE,
         std.err = TRUE, corr = FALSE, method = "Nelder-Mead")
{
    if (missing(x) || length(x) == 0 || mode(x) != "numeric") 
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
                dext(x, densfun, distnfun, p, ...)
    if(l > 1)
        body(dens) <- parse(text = paste("dext(x, densfun, distnfun,",
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
        n = length(x), model = "ext"),
        class = "evd")
}

"forder"<-
function(x, start, densfun, distnfun, ..., distn, mlen = 1, j = 1,
         largest = TRUE, std.err = TRUE, corr = FALSE, method = "Nelder-Mead")
{
    if (missing(x) || length(x) == 0 || mode(x) != "numeric") 
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
        n = length(x), model = "order"),
        class = "evd")
}

"fgumbel" <- function (...) .Defunct()
"frweibull" <- function (...) .Defunct()
"ffrechet" <- function (...) .Defunct()

"fgev"<-
function(x, start, ..., nsloc = NULL, std.err = TRUE, corr = FALSE, method = "BFGS", warn.inf = TRUE)
{
    if (missing(x) || length(x) == 0 || mode(x) != "numeric") 
        stop("`x' must be a non-empty numeric vector")
    call <- match.call()
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
            x, n,
            as.double(loc), as.double(scale),
            as.double(shape), dns = double(1),
            PACKAGE = "evd")$dns
    }
    if(!is.null(nsloc)) {
        if(is.vector(nsloc))
            nsloc <- data.frame(trend = nsloc)
        if(!is.data.frame(nsloc))
            stop("`nsloc' must be a vector or data frame")
        if(nrow(nsloc) != length(x))
            stop("`nsloc' and `x' are not compatible")
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
    structure(list(estimate = opt$par, std.err = std.err,
        fixed = unlist(fixed.param), param = param,
        deviance = 2*opt$value, corr = corr,
        convergence = opt$convergence, counts = opt$counts,
        message = opt$message,
        call = call, data = x, tdata = x2, nsloc = nsloc,
        n = length(x), model = "gev"), class = "evd")
}

"fgev.quantile"<-
function(x, start, ..., prob, nsloc = NULL, std.err = TRUE, corr = FALSE, method = "Nelder-Mead", warn.inf = TRUE)
{
    if (missing(x) || length(x) == 0 || mode(x) != "numeric") 
        stop("`x' must be a non-empty numeric vector")
    call <- match.call()
    nlgev <- function(quantile, scale, shape)
    {
        if(scale <= 0) return(1e6)
        quantile <- rep(quantile, length.out = length(x))
        if(shape == 0) loc <- quantile + scale * log(-log(1-prob))
        else loc <- quantile + scale/shape * (1 - (-log(1-prob))^(-shape))
        if(!is.null(nsloc)) {
            ns <- numeric(length(loc.param) - 1)
            for(i in 1:length(ns))
                ns[i] <- get(loc.param[i+1])
            loc <- drop(nslocmat %*% ns) + loc
        }
        .C("nlgev",
            x, n,
            as.double(loc), as.double(scale),
            as.double(shape), dns = double(1),
            PACKAGE = "evd")$dns
    }
    if(!is.null(nsloc)) {
        if(is.vector(nsloc))
            nsloc <- data.frame(trend = nsloc)
        if(!is.data.frame(nsloc))
            stop("`nsloc' must be a vector or data frame")
        if(nrow(nsloc) != length(x))
            stop("`nsloc' and `x' are not compatible")
        nsloc <- nsloc[!is.na(x), ,drop = FALSE]
        nslocmat <- as.matrix(nsloc)
    }
    x <- as.double(x[!is.na(x)])
    n <- as.integer(length(x))
    if(is.null(nsloc)) loc.param <- "quantile"
    else loc.param <- c("quantile", paste("loc", names(nsloc), sep=""))
    param <- c(loc.param, "scale", "shape")
    if(missing(start)) {
        start <- as.list(numeric(length(param)))
        names(start) <- param
        start$scale <- sqrt(6 * var(x))/pi
        start.loc <- mean(x) - 0.58 * start$scale
        start$quantile <- start.loc - start$scale * log(-log(1-prob))
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
    structure(list(estimate = opt$par, std.err = std.err,
        fixed = unlist(fixed.param), param = param,
        deviance = 2*opt$value, corr = corr, convergence = opt$convergence,
        counts = opt$counts, message = opt$message,
        call = call, data = x, tdata = x2,
        nsloc = nsloc, n = length(x), model = "gev.quantile"),
        prob = prob, loc = loc, class = "evd")
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
"deviance.evd" <- function (object, ...) object$deviance
"std.errors" <- function (object, ...) UseMethod("std.errors")
"std.errors.evd" <- function (object, ...) object$std.err

"anova.evd" <- function (object, object2, ...) 
{
    narg <- nargs()
    if(narg < 2) stop("there must be two or more arguments")
    dots <- as.list(substitute(list(...)))[-1]
    dots <- sapply(dots,function(x) deparse(x))
    if(!length(dots)) dots <- NULL
    model1 <- deparse(substitute(object))
    model2 <- deparse(substitute(object2))
    models <- c(model1, model2, dots)
    for(i in 1:narg) {
        if(!inherits(get(models[i]), "evd")) 
            stop("Use only with 'evd' objects")
    }
    for(i in 1:(narg-1)) {
        a <- names(get(models[i])$estimate)
        b <- names(get(models[i+1])$estimate)
        if(!all(b %in% a))
            stop("models are not nested")
    }
    dv <- npar <- numeric(narg)
    for(i in 1:narg) {
        dv[i] <- get(models[i])$deviance
        npar[i] <- length(get(models[i])$estimate)
    }
    df <- -diff(npar)
    if(any(df == 0)) stop("models are not nested")
    dvdiff <- diff(dv)
    if(any(dvdiff < 0)) stop("models are not nested")
    pval <- pchisq(dvdiff, df = df, lower.tail = FALSE)
    table <- data.frame(npar, dv, c(NA,df), c(NA,dvdiff), c(NA,pval))
    dimnames(table) <- list(models, c("M.Df", "Deviance", "Df", "Chisq",
                                      "Pr(>chisq)"))
    structure(table, heading = c("Analysis of Deviance Table\n"),
              class = c("anova", "data.frame"))
}

"plot.evd" <-  function(x, which = 1:4, main = c("Probability Plot",
     "Quantile Plot", "Density Plot", "Return Level Plot"),
     ask = nb.fig < length(which) && dev.interactive(), ci = TRUE,
     adjust = 1, jitter = FALSE, nplty = 2, ...) 
{
    if (!inherits(x, "evd")) 
        stop("Use only with `evd' objects")
    if (!is.numeric(which) || any(which < 1) || any(which > 4)) 
        stop("`which' must be in 1:4")
    if(x$model %in% c("ext","order"))
        stop("diagnostic plots not implemented for this model")
    show <- rep(FALSE, 4)
    show[which] <- TRUE
    nb.fig <- prod(par("mfcol"))
    if (ask) {
        op <- par(ask = TRUE)
        on.exit(par(op))
    }
    if (show[1]) {
        pp(x, ci = ci, main = main[1], xlim = c(0,1), ylim = c(0,1), ...)
    }
    if (show[2]) {
        qq(x, ci = ci, main = main[2], ...)
    }
    if (show[3]) {
        dens(x, adjust = adjust, nplty = nplty, jitter = jitter,
             main = main[3], ...)
    }
    if (show[4]) {
        rl(x, ci = ci, main = main[4], ...)
    }
    invisible(x)
}

"qq" <-  function(x, ci = TRUE, main = "Quantile Plot", xlab = "Model", ylab = "Empirical", ...)
{
    if(!(x$model %in% c("gev","gev.quantile")))
        stop("qq plot not implemented for this model")
    if(x$model == "gev") loc <- x$param["loc"]
    if(x$model == "gev.quantile") loc <- attributes(x)$loc
    quant <- qgev(ppoints(x$tdata), loc = loc,
                 scale = x$param["scale"], shape = x$param["shape"])
    if(!ci) {
      plot(quant, sort(x$tdata), main = main, xlab = xlab, ylab = ylab, ...)
      abline(0, 1)
    }
    else {
      samp <- rgev(x$n*99, loc = loc,
                 scale = x$param["scale"], shape = x$param["shape"])
      samp <- matrix(samp, x$n, 99)
      samp <- apply(samp, 2, sort)
      samp <- apply(samp, 1, sort)
      env <- t(samp[c(3,97),])
      rs <- sort(x$tdata)
      matplot(quant, cbind(rs,env), main = main, xlab = xlab, ylab = ylab,
              type = "pnn", pch = 4, ...)
      xyuser <- par("usr")
      smidge <- min(diff(c(xyuser[1], quant, xyuser[2])))/2
      segments(quant-smidge, env[,1], quant+smidge, env[,1])
      segments(quant-smidge, env[,2], quant+smidge, env[,2])
      abline(0, 1)
    }
    invisible(list(x = quant, y = sort(x$tdata)))
}

"pp" <-  function(x, ci = TRUE, main = "Probability Plot", xlab = "Empirical", ylab = "Model", ...)
{
    if(!(x$model %in% c("gev","gev.quantile")))
        stop("pp plot not implemented for this model")
    ppx <- ppoints(x$n)
    if(x$model == "gev") loc <- x$param["loc"]
    if(x$model == "gev.quantile") loc <- attributes(x)$loc
    probs <- pgev(sort(x$tdata), loc = loc,
                 scale = x$param["scale"], shape = x$param["shape"])
    if(!ci) {
        plot(ppx, probs, main = main, xlab = xlab, ylab = ylab, ...)
        abline(0, 1)
    }
    else {
        samp <- rgev(x$n*99, loc = loc,
                   scale = x$param["scale"], shape = x$param["shape"])
        samp <- matrix(samp, x$n, 99)
        samp <- apply(samp, 2, sort)
        samp <- apply(samp, 1, sort)
        env <- t(samp[c(3,97),])
        env[,1] <- pgev(env[,1], loc = loc,
                    scale = x$param["scale"], shape = x$param["shape"])
        env[,2] <- pgev(env[,2], loc = loc,
                    scale = x$param["scale"], shape = x$param["shape"])
        matplot(ppx, cbind(probs, env), main = main, xlab = xlab,
                ylab = ylab, type = "pnn", pch = 4, ...)
        xyuser <- par("usr")
        smidge <- min(diff(c(xyuser[1], ppx, xyuser[2])))/2
        segments(ppx-smidge, env[,1], ppx+smidge, env[,1])
        segments(ppx-smidge, env[,2], ppx+smidge, env[,2])
        abline(0, 1)
    }
    invisible(list(x = ppoints(x$n), y = probs))
}

"rl" <-  function(x, ci = TRUE, main = "Return Level Plot", xlab = "-1/log(1-1/Return Period)", ylab = "Return Level", ...)
{
    if(!(x$model %in% c("gev","gev.quantile")))
        stop("return level plot not implemented for this model")
    ppx <- ppoints(x$tdata)
    if(x$model == "gev") loc <- x$param["loc"]
    if(x$model == "gev.quantile") loc <- attributes(x)$loc
    rps <- c(1.001,10^(seq(0,3,len=200))[-1])
    p.upper <- 1/rps
    rlev <- qgev(p.upper, loc = loc, scale = x$param["scale"],
              shape = x$param["shape"], lower.tail = FALSE)
    if(!ci) {
        plot(-1/log(ppx), sort(x$tdata),log = "x", main = main,
             xlab = xlab, ylab = ylab, ...)
        lines(-1/log(1-p.upper), rlev)
    }
    else {
        samp <- rgev(x$n*99, loc = loc,
                   scale = x$param["scale"], shape = x$param["shape"])
        samp <- matrix(samp, x$n, 99)
        samp <- apply(samp, 2, sort)
        samp <- apply(samp, 1, sort)
        env <- t(samp[c(3,97),])
        rs <- sort(x$tdata)
        matplot(-1/log(ppx), cbind(rs,env), main = main, xlab = xlab,
                ylab = ylab, type = "pnn", pch = 4, log = "x", ...)
        lines(-1/log(1-p.upper), rlev)
        xyuser <- par("usr")
        smidge <- min(diff(c(xyuser[1], log10(-1/log(ppx)), xyuser[2])))/2
        segments((-1/log(ppx))*exp(-smidge), env[,1],
                 (-1/log(ppx))*exp(smidge), env[,1])
        segments((-1/log(ppx))*exp(-smidge), env[,2],
                 (-1/log(ppx))*exp(smidge), env[,2])
    }
    invisible(list(x = -1/log(1-p.upper), y = rlev))
}

"dens" <-  function(x, adjust = 1, nplty = 2, jitter = FALSE, main = "Density Plot", xlab = "Quantile", ylab = "Density", ...)
{
    if(!(x$model %in% c("gev","gev.quantile")))
        stop("density plot not implemented for this model")
    if(x$model == "gev") loc <- x$param["loc"]
    if(x$model == "gev.quantile") loc <- attributes(x)$loc
    xlimit <- range(x$tdata)
    xlimit[1] <- xlimit[1] - diff(xlimit) * 0.075
    xlimit[2] <- xlimit[2] + diff(xlimit) * 0.075
    xvec <- seq(xlimit[1], xlimit[2], length = 100)
    dens <- dgev(xvec, loc = loc, scale = x$param["scale"],
                shape = x$param["shape"])
    plot(spline(xvec, dens), main = main, xlab = xlab, ylab = ylab,
         type = "l", ...)
    if(jitter) rug(jitter(x$tdata))
    else rug(x$tdata)
    lines(density(x$tdata, adjust = adjust), lty = nplty)
    invisible(list(x = xvec, y = dens))
}

"profile.evd" <-  function(fitted, which = names(fitted$estimate), conf = 0.999, mesh = fitted$std.err[which]/2, xmin = rep(-Inf, length(which)), xmax = rep(Inf, length(which)), convergence = FALSE, control = list(maxit = 5000), ...)
{
    if (!inherits(fitted, "evd")) 
        stop("Use only with `evd' objects")
    if(fitted$model %in% c("ext","order"))
        stop("profiles not implemented for this model")
    if(length(xmin) != length(which))
        stop("`xmin' and `which' must have the same length")
    if(length(xmax) != length(which))
        stop("`xmax' and `which' must have the same length")
    if(length(fitted$estimate) < 2)
        stop("cannot profile one dimensional likelihood")
    if(!is.character(which))
        stop("`which' must be a character vector")
    if(!all(which %in% names(fitted$estimate)))
        stop("`which' contains unrecognized or unestimated parameters")
    if(is.null(fitted$std.err) && missing(mesh))
       stop("fitted model must contain standard errors")
    prof.list <- as.list(numeric(length(which)))
    names(xmin) <- names(xmax) <- names(prof.list) <- which
    mles <- fitted$estimate[which]                   
    for(j in which) {
        print(paste("profiling",j))
        prof <- matrix(NA, nrow = 64, ncol = length(fitted$estimate) + 1)
        parvec1 <- seq(mles[j] + mesh[j], mles[j] + 16*mesh[j], length = 32) 
        parvec2 <- seq(mles[j] - mesh[j], mles[j] - 16*mesh[j], length = 32)
        if(any(parvec1 >= xmax[j]))
            parvec1 <- c(parvec1[parvec1 < xmax[j]], xmax[j])
        if(any(parvec2 <= xmin[j]))
            parvec2 <- c(parvec2[parvec2 > xmin[j]], xmin[j])
        start <- as.list(fitted$estimate[!names(fitted$estimate) %in% j])
        call.args <- c(list(fitted$data, start, 0), as.list(fitted$fixed),
           list(FALSE, FALSE, "Nelder-Mead", FALSE, control))
        names(call.args) <- c("x", "start", j, names(fitted$fixed),
           "std.err", "corr", "method", "warn.inf", "control")
        dimnames(prof) <- list(NULL, c(j, "deviance", names(start)))
        call.fn <- paste("f", fitted$model, sep="")
        if(!is.na(pmatch("gev", fitted$model)))
            call.args$nsloc <- fitted$nsloc
        if(!is.na(pmatch("bv", fitted$model))) {
            call.args$nsloc1 <- fitted$nsloc1
            call.args$nsloc2 <- fitted$nsloc2
        }
        if(fitted$model == "gev.quantile")
            call.args$prob <- attributes(fitted)$prob
        for(i in 1:32) {
            call.args[[j]] <- parvec1[i]
            fit.mod <- do.call(call.fn, call.args)
            if(convergence) print(fit.mod$convergence)
            call.args[["start"]] <- as.list(fit.mod$estimate)
            prof[i+32,1] <- parvec1[i]
            prof[i+32,2] <- fit.mod$deviance
            prof[i+32,-(1:2)] <- fit.mod$estimate
            if(abs(fit.mod$deviance - fitted$deviance) > qchisq(0.999,1))
              break;
            if(parvec1[i] == xmax[j]) break;
        }
        for(i in 1:32) {
            call.args[[j]] <- parvec2[i]
            call.args[["start"]] <-
              as.list(fitted$estimate[!names(fitted$estimate) %in% j])
            fit.mod <- do.call(call.fn, call.args)
            if(convergence) print(fit.mod$convergence)
            call.args[["start"]] <- as.list(fit.mod$estimate)
            prof[33-i,1] <- parvec2[i]
            prof[33-i,2] <- fit.mod$deviance
            prof[33-i,-(1:2)] <- fit.mod$estimate
            if(abs(fit.mod$deviance - fitted$deviance) > qchisq(0.999,1))
              break;
            if(parvec2[i] == xmin[j]) break;
        }
        prof <- na.omit(prof)
        attributes(prof)$na.action <- NULL
        prof.list[[j]] <- prof
    }
    structure(prof.list, deviance = fitted$deviance,
              xmin = xmin, xmax = xmax, class = "profile.evd")
}

"pcint" <- function(prof, which = names(prof), ci = 0.95)
{
    if (!inherits(prof, "profile.evd")) 
        stop("Use only with `profile.evd' objects")
    if(!is.character(which))
        stop("`which' must be a character vector")
    if(!all(which %in% names(prof)))
        stop("`which' contains unprofiled parameters")
    rdev <- attributes(prof)$deviance + qchisq(ci, df = 1)
    ci <- as.list(which)
    names(ci) <- which
    for(i in which) {
        x <- prof[[i]]
        n <- nrow(x)
        ulim <- min(x[1,"deviance"],x[n,"deviance"])
        llim <- min(x[1,"deviance"],x[n,"deviance"])
        th.l <- x[1, i] == attributes(prof)$xmin[i]
        th.u <- x[n, i] == attributes(prof)$xmax[i]
        if(x[1,"deviance"] <= rdev && !th.l)
           stop("confidence coefficient is too high")
        if(x[n,"deviance"] <= rdev && !th.u)
           stop("confidence coefficient is too high")
        halves <- c(diff(x[,"deviance"]) < 0, FALSE)
        if(x[1,"deviance"] <= rdev && th.l)
            lower <- x[1, i]
        else lower <- approx(x[halves,2], x[halves,1], xout = rdev)$y
        if(x[n,"deviance"] <= rdev && th.u)
            upper <- x[n, i]
        else upper <- approx(x[!halves,2], x[!halves,1], xout = rdev)$y
        ci[[i]] <- c(lower,upper)
    }
    ci
}
    
profile2d <- function (fitted, ...) {
UseMethod("profile2d")
}

"profile2d.evd" <-  function(fitted, prof, which, pts = 20, convergence = FALSE, control = list(maxit = 5000), ...)
{
    if(fitted$model %in% c("ext","order"))
        stop("profiles not implemented for this model")
    if (!inherits(prof, "profile.evd")) 
        stop("`prof' must be a `profile.evd' object")
    if(length(fitted$estimate) < 3)
        stop("Cannot profile two dimensional likelihood")
    if(missing(which) || !is.character(which) || length(which) != 2)
        stop("`which' must be a character vector of length two")
    if(!all(which %in% names(fitted$estimate)))
        stop("`which' contains unrecognized or unestimated parameters")
    if(!all(which %in% names(prof)))
        stop("`which' contains unprofiled parameters")
    if(is.null(fitted$std.err))
       stop("fitted model must contain standard errors")
    prof.list <- as.list(numeric(3))
    names(prof.list) <- c("trace",which)
    limits1 <- range(prof[[which[1]]][,1])
    limits2 <- range(prof[[which[2]]][,1])
    mles <- fitted$estimate[which]                   
    prof <- matrix(NA, nrow = pts^2, ncol = length(fitted$estimate) + 1)
    parvec1 <- seq(limits1[1], limits1[2], length = pts)
    prof.list[[which[1]]] <- parvec1
    parvec2 <- seq(limits2[1], limits2[2], length = pts)
    prof.list[[which[2]]] <- parvec2
    pars <- expand.grid(parvec1, parvec2)
    start <- as.list(fitted$estimate[!names(fitted$estimate) %in% which])
    call.args <- c(list(fitted$data, start, 0, 0), as.list(fitted$fixed),
       list(FALSE, FALSE, "Nelder-Mead", FALSE, control))
    names(call.args) <- c("x", "start", which[1], which[2],
       names(fitted$fixed), "std.err", "corr", "method",
       "warn.inf", "control")
    dimnames(prof) <- list(NULL, c(which, "deviance", names(start)))
    call.fn <- paste("f", fitted$model, sep="")
    if(!is.na(pmatch("gev", fitted$model)))
            call.args$nsloc <- fitted$nsloc
    if(!is.na(pmatch("bv", fitted$model))) {
            call.args$nsloc1 <- fitted$nsloc1
            call.args$nsloc2 <- fitted$nsloc2
    }
    if(fitted$model == "gev.quantile")
        call.args$prob <- attributes(fitted)$prob
    for(i in 1:pts^2) {
        call.args[[which[1]]] <- pars[i,1]
        call.args[[which[2]]] <- pars[i,2]
        fit.mod <- do.call(call.fn, call.args)
        if(convergence) print(fit.mod$convergence)
        prof[i,1] <- pars[i,1]
        prof[i,2] <- pars[i,2]
        prof[i,3] <- fit.mod$deviance
        prof[i,-(1:3)] <- fit.mod$estimate
    }
    prof.list[["trace"]] <- prof
    if(any(prof[,"deviance"] == 2e6))
        warning("non-convergence present in profile2d object")
    structure(prof.list, deviance = fitted$deviance, class = "profile2d.evd")
}

"plot.profile.evd" <-  function(x, which = names(x), main = NULL,
     ask = nb.fig < length(which) && dev.interactive(), ci = 0.95,
     clty = 2, ...) 
{
    if (!inherits(x, "profile.evd")) 
        stop("Use only with `profile.evd' objects")
    if(!is.character(which))
        stop("`which' must be a character vector")
    if(!all(which %in% names(x)))
        stop("`which' contains unprofiled parameters")
    nb.fig <- prod(par("mfcol"))
    if (ask) {
        op <- par(ask = TRUE)
        on.exit(par(op))
    }
    if(is.null(main))
        main <- paste("Profile Deviance of", which, "Parameter")
    for(i in which) {
        plot(spline(x[[i]][,1], x[[i]][,2], n = 75), type = "l",
             xlab = i, ylab = "profile deviance",
             main = main[match(i,which)], ...)
        cdist <- attributes(x)$deviance + qchisq(ci, df = 1)
        abline(h = cdist, lty = clty)
    }
    invisible(x)
}

"plot.profile2d.evd" <-  function(x, main = NULL, ci = c(0.5,0.8,0.9,0.95,0.975, 0.99, 0.995), col = heat.colors(8), intpts = 75, ...) 
{
    if (!inherits(x, "profile2d.evd")) 
        stop("Use only with `profile2d.evd' objects")
    which <- names(x)[2:3]
    if(is.null(main))
        main <- paste("Profile Deviance of", which[1],
                      "and", which[2])
    br.pts <- attributes(x)$deviance + qchisq(c(0,ci), df = 2)
    prof <- x$trace[,"deviance"]
    if(any(prof == 2e6))
        warning("non-convergence present in profile2d object")
    if(!require(akima)) {
        image(x[[which[1]]], x[[which[2]]],
              matrix(prof, nrow = length(x[[which[1]]])),
              col = col, breaks = c(br.pts, 2e6+1),
              main = main, xlab = which[1], ylab = which[2], ...)
    }
    else {
        lim1 <- range(x[[which[1]]])
        lim2 <- range(x[[which[2]]]) 
        prof.interp <- interp(x$trace[,1], x$trace[,2], prof,
            xo = seq(lim1[1], lim1[2], length = intpts),
            yo = seq(lim2[1], lim2[2], length = intpts))
        image(prof.interp, col = col, breaks = c(br.pts, max(prof)),
              main = main, xlab = which[1], ylab = which[2], ...)
    }
    invisible(x)
}







