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
    if(min(p) <= 0 || max(p) >=1)
        stop("`p' must contain probabilities in (0,1)")
    if(min(scale) < 0 || min(shape) <= 0) stop("invalid arguments")
    if(!lower.tail) p <- 1 - p
    loc - scale * (-log(p))^(1/shape)
}

"qgev"<-
function(p, loc = 0, scale = 1, shape = 0, lower.tail = TRUE)
{
    if(min(p) <= 0 || max(p) >=1)
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
    if(min(p) <= 0 || max(p) >=1)
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
    xpos <- x[x>0]
    nn <- length(x)
    scale <- rep(scale, length.out = nn)[x>0]
    shape <- rep(shape, length.out = nn)[x>0]
    d <- numeric(nn)
    d[x>0] <- log(shape/scale) - (1+shape) * log(xpos) - xpos^(-shape)
    d[x<=0] <- -Inf
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
    xneg <- x[x<0]
    nn <- length(x)
    scale <- rep(scale, length.out = nn)[x<0]
    shape <- rep(shape, length.out = nn)[x<0]
    d <- numeric(nn)
    d[x<0] <- log(shape/scale) + (shape-1) * log(-xneg) - (-xneg)^shape
    d[x>=0] <- -Inf
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
        xxpos <- xx[xx>0]
        scale <- rep(scale, length.out = nn)[xx>0]
        d <- numeric(nn)
        d[xx>0] <- log(1/scale) - xxpos^(-1/shape) - (1/shape + 1)*log(xxpos)
        d[xx<=0] <- -Inf
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

"fdensfun"<-
function(x, densfun, start, ...)
{
    if (missing(x) || length(x) == 0) 
        stop("`x' must be a non-empty object")
    if (is.null(start) || !is.list(start)) 
        stop("`start' must be a named list")
    nllh <- function(p, ...) {
        dvec <- dens(p, ..., log = TRUE)
        if(any(is.infinite(dvec)))
            return(1e6)
        else 
            return(-sum(dvec))
    }
    nm <- names(start)
    l <- length(nm)
    f <- formals(densfun)
    args <- names(f)
    m <- match(nm, args)
    if (any(is.na(m))) 
        stop("`start' specifies unknown arguments")
    formals(densfun) <- c(f[c(1, m)], f[-c(1, m)])
    dens <- function(p, x, ...) densfun(x, p, ...)
    if(l > 1)
        body(dens) <- parse(text = paste("densfun(x,", paste("p[",1:l,
            "]", collapse = ", "), ", ...)"))
    fixed.par <- list(...)
    fixed.par <- fixed.par[names(fixed.par) %in% args]
    start.arg <- c(list(p = unlist(start), x = x), fixed.par)
    if(do.call("nllh", start.arg) == 1e6)
        warning("negative log-likelihood is infinite at starting values")
    opt <- optim(start, nllh, x = x, hessian = TRUE, ...)
    if (opt$convergence > 0) 
        warning("optimization may not have succeeded")
    opt$std.err <- sqrt(diag(solve(opt$hessian)))
    #opt$std.err <- rep(NA,l)
    names(opt$std.err) <- nm
    list(estimate = opt$par, std.err = opt$std.err, deviance = 2*opt$value,
         counts = opt$counts)
}

"fext"<-
function(x, start, densfun, distnfun, ..., distn, mlen = 1, largest = TRUE)
{
    if (missing(x) || length(x) == 0 || mode(x) != "numeric") 
        stop("`x' must be a non-empty numeric object")
    if (is.null(start) || !is.list(start)) 
        stop("`start' must be a named list")
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
                 largest = largest)
    if (opt$convergence > 0) 
        stop("optimization failed")
    opt$std.err <- sqrt(diag(solve(opt$hessian)))
    #opt$std.err <- rep(NA,l)
    names(opt$std.err) <- nm
    list(estimate = opt$par, std.err = opt$std.err, deviance = 2*opt$value,
         counts = opt$counts)
}

"forder"<-
function(x, start, densfun, distnfun, ..., distn, mlen = 1, j = 1,
         largest = TRUE)
{
    if (missing(x) || length(x) == 0 || mode(x) != "numeric") 
        stop("`x' must be a non-empty numeric object")
    if (is.null(start) || !is.list(start)) 
        stop("`start' must be a named list")
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
                 distnfun = distnfun, mlen = mlen, j = j, largest = largest)
    if (opt$convergence > 0) 
        stop("optimization failed")
    opt$std.err <- sqrt(diag(solve(opt$hessian)))
    #opt$std.err <- rep(NA,l)
    names(opt$std.err) <- nm
    list(estimate = opt$par, std.err = opt$std.err, deviance = 2*opt$value,
         counts = opt$counts)
}

"ffrechet"<-
function(x, start, ...)
{
    ddfrechet <- function(x, loc = 0, scale = 1, shape = 1, log) {
                     if(any(c(scale, shape) <= 0))
                         return(-Inf)
                     dfrechet(x = x, loc = loc, scale = scale,
                              shape = shape, log = log)
                 }
    fdensfun(x = x, densfun = ddfrechet, start = start, ...)
}

"fgumbel"<-
function(x, start, ...)
{
    fgev(x = x, start = start, shape = 0, ...)
}

"frweibull"<-
function(x, start, ...)
{
    ddrweibull <- function(x, loc = 0, scale = 1, shape = 1, log) {
                      if(any(c(scale, shape) <= 0))
                          return(-Inf)
                      drweibull(x = x, loc = loc, scale = scale,
                                shape = shape, log = log)
                  }    
    fdensfun(x = x, densfun = ddrweibull, start = start, ...)
}

"fgev"<-
function(x, start, ...)
{
    ddgev <- function(x, loc = 0, scale = 1, shape = 0, log) {
                 if(scale <= 0) return(-Inf)
                 dgev(x = x, loc = loc, scale = scale,
                      shape = shape, log = log)
                 } 
    fdensfun(x = x, densfun = ddgev, start = start, ...)
}











