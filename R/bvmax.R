"rbvlog"<-
# Uses Algorithm 1.1 in Stephenson(2003)
function(n, dep, mar1 = c(0,1,0), mar2 = mar1)
{
    if(length(dep) != 1 || mode(dep) != "numeric" || dep <= 0 ||
        dep > 1) stop("invalid argument for `dep'")
    sim <- .C("rbvlog_shi",
               as.integer(n), as.double(dep), sim = double(2*n),
               PACKAGE = "evd")$sim
    sim <- matrix(sim, nrow = n, ncol = 2, byrow = TRUE)
    mtransform(1/sim, list(mar1, mar2), inv = TRUE, drp = TRUE)
}

"rbvalog"<-
# Uses Algorithm 1.2 in Stephenson(2003)
function(n, dep, asy = c(1,1), mar1 = c(0,1,0), mar2 = mar1)
{
    if(length(dep) != 1 || mode(dep) != "numeric" || dep <= 0 ||
       dep > 1) stop("invalid argument for `dep'")
    if(length(asy) != 2 || mode(asy) != "numeric" || min(asy) < 0 ||
       max(asy) > 1) stop("invalid argument for `asy'")
    if(dep == 1 || any(asy == 0)) {
        asy <- c(0,0)
        dep <- 1
    }
    sim <- .C("rbvalog_shi",
              as.integer(n), as.double(dep), as.double(asy),
              sim = double(2*n), PACKAGE = "evd")$sim
    sim <- matrix(sim, nrow = n, ncol = 2, byrow = TRUE)
    mtransform(1/sim, list(mar1, mar2), inv = TRUE, drp = TRUE)
}

"rbvhr" <-
function(n, dep, mar1 = c(0,1,0), mar2 = mar1)
{
    if(length(dep) != 1 || mode(dep) != "numeric" || dep <= 0)
        stop("invalid argument for `dep'")
    sim <- .C("rbvhr",
        as.integer(n), as.double(dep), sim = runif(2*n),
        PACKAGE = "evd")$sim
    sim <- matrix(sim, nrow = n, ncol = 2, byrow = TRUE)
    mtransform(-log(sim), list(mar1, mar2), inv = TRUE, drp = TRUE)
}

"rbvneglog"<- 
function(n, dep, mar1 = c(0,1,0), mar2 = mar1)
{
    if(length(dep) != 1 || mode(dep) != "numeric" || dep <= 0)
        stop("invalid argument for `dep'")
    sim <- .C("rbvneglog",
        as.integer(n), as.double(dep), sim = runif(2*n),
        PACKAGE = "evd")$sim
    sim <- matrix(sim, nrow = n, ncol = 2, byrow = TRUE)
    mtransform(-log(sim), list(mar1, mar2), inv = TRUE, drp = TRUE)
}

"rbvaneglog"<- 
function(n, dep, asy = c(1,1), mar1 = c(0,1,0), mar2 = mar1)
{
    if(length(dep) != 1 || mode(dep) != "numeric" || dep <= 0)
        stop("invalid argument for `dep'")
    if(length(asy) != 2 || mode(asy) != "numeric" || min(asy) < 0 ||
       max(asy) > 1) stop("invalid argument for `asy'")
    sim <- .C("rbvaneglog",
        as.integer(n), as.double(dep), as.double(asy), sim = runif(2*n),
        PACKAGE = "evd")$sim
    sim <- matrix(sim, nrow = n, ncol = 2, byrow = TRUE)
    mtransform(-log(sim), list(mar1, mar2), inv = TRUE, drp = TRUE)
}

"rbvbilog"<- 
function(n, alpha, beta, mar1 = c(0,1,0), mar2 = mar1)
{
    if(length(alpha) != 1 || mode(alpha) != "numeric")
        stop("invalid argument for `alpha'")
    if(length(beta) != 1 || mode(beta) != "numeric")
        stop("invalid argument for `beta'")
    if(any(c(alpha,beta) <= 0) || any(c(alpha,beta) >= 1))
        stop("`alpha' and `beta' must be in the open interval (0,1)")
    sim <- .C("rbvbilog",
        as.integer(n), as.double(alpha), as.double(beta), sim = runif(2*n),
        PACKAGE = "evd")$sim
    sim <- matrix(sim, nrow = n, ncol = 2, byrow = TRUE)
    mtransform(-log(sim), list(mar1, mar2), inv = TRUE, drp = TRUE)
}

"rbvnegbilog"<- 
function(n, alpha, beta, mar1 = c(0,1,0), mar2 = mar1)
{
    if(length(alpha) != 1 || mode(alpha) != "numeric")
        stop("invalid argument for `alpha'")
    if(length(beta) != 1 || mode(beta) != "numeric")
        stop("invalid argument for `beta'")
    if(any(c(alpha,beta) <= 0))
        stop("`alpha' and `beta' must be non-negative")
    sim <- .C("rbvnegbilog",
        as.integer(n), as.double(alpha), as.double(beta), sim = runif(2*n),
        PACKAGE = "evd")$sim
    sim <- matrix(sim, nrow = n, ncol = 2, byrow = TRUE)
    mtransform(-log(sim), list(mar1, mar2), inv = TRUE, drp = TRUE)
}

"rbvct" <-
function(n, alpha, beta, mar1 = c(0,1,0), mar2 = mar1)
{
    if(length(alpha) != 1 || mode(alpha) != "numeric")
        stop("invalid argument for `alpha'")
    if(length(beta) != 1 || mode(beta) != "numeric")
        stop("invalid argument for `beta'")
    if(any(c(alpha,beta) <= 0))
        stop("`alpha' and `beta' must be non-negative")
    sim <- .C("rbvct",
        as.integer(n), as.double(alpha), as.double(beta), sim = runif(2*n),
        PACKAGE = "evd")$sim
    sim <- matrix(sim, nrow = n, ncol = 2, byrow = TRUE)
    mtransform(-log(sim), list(mar1, mar2), inv = TRUE, drp = TRUE)
}

"rbvevd" <-
function(n, dep, asy = c(1,1), alpha, beta, model = c("log", "alog",
    "hr", "neglog", "aneglog", "bilog", "negbilog", "ct"),
    mar1 = c(0,1,0), mar2 = mar1)
{
  model <- match.arg(model)
  m1 <- c("bilog", "negbilog", "ct")
  m2 <- c(m1, "log", "hr", "neglog")
  m3 <- c("log", "alog", "hr", "neglog", "aneglog")
  if((model %in% m1) && !missing(dep))
    warning("ignoring `dep' argument")
  if((model %in% m2) && !missing(asy))
    warning("ignoring `asy' argument")
  if((model %in% m3) && !missing(alpha))
    warning("ignoring `alpha' argument")
  if((model %in% m3) && !missing(beta))
    warning("ignoring `beta' argument")
    
  switch(model,
    log = rbvlog(n = n, dep = dep, mar1 = mar1, mar2 = mar2),
    alog = rbvalog(n = n, dep = dep, asy = asy, mar1 = mar1, mar2 = mar2),
    hr = rbvhr(n = n, dep = dep, mar1 = mar1, mar2 = mar2),
    neglog = rbvneglog(n = n, dep = dep, mar1 = mar1, mar2 = mar2),
    aneglog = rbvaneglog(n = n, dep = dep, asy = asy, mar1 = mar1,
      mar2 = mar2),
    bilog = rbvbilog(n = n, alpha = alpha, beta = beta, mar1 = mar1,
      mar2 = mar2),
    negbilog = rbvnegbilog(n = n, alpha = alpha, beta = beta, mar1 = mar1,
      mar2 = mar2),
    ct = rbvct(n = n, alpha = alpha, beta = beta, mar1 = mar1, mar2 = mar2)) 
}

"evmc" <-
function(n, dep, asy = c(1,1), alpha, beta, model = c("log", "alog",
    "hr", "neglog", "aneglog", "bilog", "negbilog", "ct"),
    margins = c("uniform","exponential","frechet","gumbel"))
{
  model <- match.arg(model)
  m1 <- c("bilog", "negbilog", "ct")
  m2 <- c(m1, c("log", "hr", "neglog"))
  m3 <- c("log", "alog", "hr", "neglog", "aneglog")
  if((model %in% m1) && !missing(dep))
    warning("ignoring `dep' argument")
  if((model %in% m2) && !missing(asy))
    warning("ignoring `asy' argument")
  if((model %in% m3) && !missing(alpha))
    warning("ignoring `alpha' argument")
  if((model %in% m3) && !missing(beta))
    warning("ignoring `beta' argument")

  nn <- as.integer(1)
  if(!(model %in% m1)) {
    if(length(dep) != 1 || mode(dep) != "numeric" || dep <= 0)
        stop("invalid argument for `dep'")
    if((model %in% c("log", "alog")) && dep > 1)
        stop("`dep' must be in the interval (0,1]")
    dep <- as.double(dep)
  }
  if(!(model %in% m2)) {
    if(length(asy) != 2 || mode(asy) != "numeric" || min(asy) < 0 ||
       max(asy) > 1) stop("invalid argument for `asy'")
    if(model == "alog" && (dep == 1 || any(asy == 0))) {
        asy <- c(0,0)
        dep <- 1
    }
    asy <- as.double(asy[c(2,1)])
  }
  if(!(model %in% m3)) {
    if(length(alpha) != 1 || mode(alpha) != "numeric" || alpha <= 0)
        stop("invalid argument for `alpha'")
    if(length(beta) != 1 || mode(beta) != "numeric" || beta <= 0)
        stop("invalid argument for `beta'")
    if(model == "bilog" && any(c(alpha,beta) >= 1))
        stop("`alpha' and `beta' must be in the open interval (0,1)")
    alpha <- as.double(beta)
    beta <- as.double(alpha)
  } 
    
  evmc <- runif(n)
  for(i in 2:n) {
    evmc[c(i,i-1)] <- switch(model,
      log = .C("rbvlog", nn, dep, sim = evmc[c(i,i-1)], PACKAGE = "evd")$sim,
      alog = .C("rbvalog", nn, dep, asy, sim = evmc[c(i,i-1)],
        PACKAGE = "evd")$sim,
      hr = .C("rbvhr", nn, dep, sim = evmc[c(i,i-1)], PACKAGE = "evd")$sim,
      neglog = .C("rbvneglog", nn, dep, sim = evmc[c(i,i-1)],
        PACKAGE = "evd")$sim,
      aneglog = .C("rbvaneglog", nn, dep, asy, sim = evmc[c(i,i-1)],
        PACKAGE = "evd")$sim,
      bilog = .C("rbvbilog", nn, alpha, beta, sim = evmc[c(i,i-1)],
        PACKAGE = "evd")$sim,
      negbilog = .C("rbvnegbilog", nn, alpha, beta, sim = evmc[c(i,i-1)],
        PACKAGE = "evd")$sim,
      ct = .C("rbvct", nn, alpha, beta, sim = evmc[c(i,i-1)],
        PACKAGE = "evd")$sim)
  }

  switch(match.arg(margins),
    frechet = -1/log(evmc), uniform = evmc,
    exponential = -log(evmc), gumbel = -log(-log(evmc))) 
}

"pbvlog"<- 
function(q, dep, mar1 = c(0,1,0), mar2 = mar1, lower.tail = TRUE)
{
    if(length(dep) != 1 || mode(dep) != "numeric" || dep <= 0 ||
        dep > 1) stop("invalid argument for `dep'")
    if(is.null(dim(q))) dim(q) <- c(1,2)
    q <- mtransform(q, list(mar1, mar2))
    v <- apply(q^(1/dep),1,sum)^dep
    pp <- exp(-v)
    if(!lower.tail) {
      pp <- 1 - pgev(-log(q[,1])) - pgev(-log(q[,2])) + pp
    }
    pp
}

"pbvalog"<- 
function(q, dep, asy = c(1,1), mar1 = c(0,1,0), mar2 = mar1, lower.tail = TRUE)
{
    if(length(dep) != 1 || mode(dep) != "numeric" || dep <= 0 ||
        dep > 1) stop("invalid argument for `dep'")
    if(length(asy) != 2 || mode(asy) != "numeric" || min(asy) < 0 ||
       max(asy) > 1) stop("invalid argument for `asy'")
    if(is.null(dim(q))) dim(q) <- c(1,2)
    q <- mtransform(q, list(mar1, mar2))
    asy <- rep(asy,rep(nrow(q),2))
    v <- apply((asy*q)^(1/dep),1,sum)^dep + apply((1-asy)*q,1,sum)
    pp <- exp(-v)
    if(!lower.tail) {
      pp <- 1 - pgev(-log(q[,1])) - pgev(-log(q[,2])) + pp
    }
    pp
}

"pbvhr" <-
function(q, dep, mar1 = c(0,1,0), mar2 = mar1, lower.tail = TRUE)
{
    if(length(dep) != 1 || mode(dep) != "numeric" || dep <= 0)
        stop("invalid argument for `dep'")
    if(is.null(dim(q))) dim(q) <- c(1,2)
    q <- mtransform(q, list(mar1, mar2))
    fn <- function(x1,x2) x1*pnorm(1/dep + dep * log(x1/x2) / 2)
    v <- fn(q[,1],q[,2]) + fn(q[,2],q[,1])
    pp <- exp(-v)
    if(!lower.tail) {
      pp <- 1 - pgev(-log(q[,1])) - pgev(-log(q[,2])) + pp
    }
    pp
}

"pbvneglog"<- 
function(q, dep, mar1 = c(0,1,0), mar2 = mar1, lower.tail = TRUE)
{
    if(length(dep) != 1 || mode(dep) != "numeric" || dep <= 0)
        stop("invalid argument for `dep'")
    if(is.null(dim(q))) dim(q) <- c(1,2)
    q <- mtransform(q, list(mar1, mar2))
    v <- apply(q,1,sum) - apply(q^(-dep),1,sum)^(-1/dep)
    pp <- exp(-v)
    if(!lower.tail) {
      pp <- 1 - pgev(-log(q[,1])) - pgev(-log(q[,2])) + pp
    }
    pp
}

"pbvaneglog"<- 
function(q, dep, asy = c(1,1), mar1 = c(0,1,0), mar2 = mar1, lower.tail = TRUE)
{
    if(length(dep) != 1 || mode(dep) != "numeric" || dep <= 0)
        stop("invalid argument for `dep'")
    if(length(asy) != 2 || mode(asy) != "numeric" || min(asy) < 0 ||
       max(asy) > 1) stop("invalid argument for `asy'")
    if(is.null(dim(q))) dim(q) <- c(1,2)
    q <- mtransform(q, list(mar1, mar2))
    asy <- rep(asy,rep(nrow(q),2))
    v <- apply(q,1,sum) - apply((asy*q)^(-dep),1,sum)^(-1/dep)
    pp <- exp(-v)
    if(!lower.tail) {
      pp <- 1 - pgev(-log(q[,1])) - pgev(-log(q[,2])) + pp
    }
    pp
}

"pbvbilog"<- 
function(q, alpha, beta, mar1 = c(0,1,0), mar2 = mar1, lower.tail = TRUE)
{
    if(length(alpha) != 1 || mode(alpha) != "numeric")
        stop("invalid argument for `alpha'")
    if(length(beta) != 1 || mode(beta) != "numeric")
        stop("invalid argument for `beta'")
    if(any(c(alpha,beta) <= 0) || any(c(alpha,beta) >= 1))
        stop("`alpha' and `beta' must be in the open interval (0,1)")
    if(is.null(dim(q))) dim(q) <- c(1,2)
    q <- mtransform(q, list(mar1, mar2))
    gma <- numeric(nrow(q))
    for(i in 1:nrow(q)) {
        gmafn <- function(x)
            (1-alpha) * q[i,1] * (1-x)^beta - (1-beta) * q[i,2] * x^alpha
        if(any(is.na(q[i,]))) gma[i] <- NA
        else if(any(is.infinite(q[i,]))) gma[i] <- 0.5
        else if(q[i,1] == 0) gma[i] <- 0
        else if(q[i,2] == 0) gma[i] <- 1
        else gma[i] <- uniroot(gmafn, lower = 0, upper = 1,
                               tol = .Machine$double.eps^0.5)$root
    }
    v <- q[,1] * gma^(1-alpha) + q[,2] * (1 - gma)^(1-beta)
    pp <- exp(-v)
    if(!lower.tail) {
      pp <- 1 - pgev(-log(q[,1])) - pgev(-log(q[,2])) + pp
    }
    pp
}

"pbvnegbilog"<- 
function(q, alpha, beta, mar1 = c(0,1,0), mar2 = mar1, lower.tail = TRUE)
{
    if(length(alpha) != 1 || mode(alpha) != "numeric")
        stop("invalid argument for `alpha'")
    if(length(beta) != 1 || mode(beta) != "numeric")
        stop("invalid argument for `beta'")
    if(any(c(alpha,beta) <= 0))
        stop("`alpha' and `beta' must be non-negative")
    if(is.null(dim(q))) dim(q) <- c(1,2)
    q <- mtransform(q, list(mar1, mar2))
    gma <- numeric(nrow(q))
    for(i in 1:nrow(q)) {
        gmafn <- function(x)
            (1+alpha) * q[i,1] * x^alpha - (1+beta) * q[i,2] * (1-x)^beta
        if(any(is.na(q[i,]))) gma[i] <- NA
        else if(any(is.infinite(q[i,]))) gma[i] <- Inf
        else if(q[i,1] == 0) gma[i] <- 1
        else if(q[i,2] == 0) gma[i] <- 0
        else gma[i] <- uniroot(gmafn, lower = 0, upper = 1,
                               tol = .Machine$double.eps^0.5)$root
    }
    v <- q[,1] + q[,2] - q[,1] * gma^(1+alpha) - q[,2] * (1 - gma)^(1+beta)
    v[is.infinite(gma)] <- Inf
    pp <- exp(-v)
    if(!lower.tail) {
      pp <- 1 - pgev(-log(q[,1])) - pgev(-log(q[,2])) + pp
    }
    pp
}

"pbvct" <-
function(q, alpha, beta, mar1 = c(0,1,0), mar2 = mar1, lower.tail = TRUE)
{
    if(length(alpha) != 1 || mode(alpha) != "numeric")
        stop("invalid argument for `alpha'")
    if(length(beta) != 1 || mode(beta) != "numeric")
        stop("invalid argument for `beta'")
    if(any(c(alpha,beta) <= 0))
        stop("`alpha' and `beta' must be non-negative")
    if(is.null(dim(q))) dim(q) <- c(1,2)
    q <- mtransform(q, list(mar1, mar2))
    u <- (alpha * q[,2]) / (alpha * q[,2] + beta * q[,1])  
    v <- q[,2] * pbeta(u, shape1 = alpha, shape2 = beta + 1) +
      q[,1] * pbeta(u, shape1 = alpha + 1, shape2 = beta, lower.tail = FALSE)
    v[is.infinite(q[,1]) || is.infinite(q[,2])] <- Inf
    v[(q[,1] + q[,2]) == 0] <- 0
    pp <- exp(-v)
    if(!lower.tail) {
      pp <- 1 - pgev(-log(q[,1])) - pgev(-log(q[,2])) + pp
    }
    pp
}

"pbvevd" <-
function(q, dep, asy = c(1,1), alpha, beta, model = c("log", "alog",
    "hr", "neglog", "aneglog", "bilog", "negbilog", "ct"),
    mar1 = c(0,1,0), mar2 = mar1, lower.tail = TRUE)
{
  model <- match.arg(model)
  m1 <- c("bilog", "negbilog", "ct")
  m2 <- c(m1, "log", "hr", "neglog")
  m3 <- c("log", "alog", "hr", "neglog", "aneglog")
  if((model %in% m1) && !missing(dep))
    warning("ignoring `dep' argument")
  if((model %in% m2) && !missing(asy))
    warning("ignoring `asy' argument")
  if((model %in% m3) && !missing(alpha))
    warning("ignoring `alpha' argument")
  if((model %in% m3) && !missing(beta))
    warning("ignoring `beta' argument")
    
  switch(model,
    log = pbvlog(q = q, dep = dep, mar1 = mar1, mar2 = mar2,
      lower.tail = lower.tail),
    alog = pbvalog(q = q, dep = dep, asy = asy, mar1 = mar1, mar2 = mar2,
      lower.tail = lower.tail),
    hr = pbvhr(q = q, dep = dep, mar1 = mar1, mar2 = mar2,
      lower.tail = lower.tail),
    neglog = pbvneglog(q = q, dep = dep, mar1 = mar1, mar2 = mar2,
      lower.tail = lower.tail),
    aneglog = pbvaneglog(q = q, dep = dep, asy = asy, mar1 = mar1,
      mar2 = mar2, lower.tail = lower.tail),
    bilog = pbvbilog(q = q, alpha = alpha, beta = beta, mar1 = mar1,
      mar2 = mar2, lower.tail = lower.tail),
    negbilog = pbvnegbilog(q = q, alpha = alpha, beta = beta, mar1 = mar1,
      mar2 = mar2, lower.tail = lower.tail),
    ct = pbvct(q = q, alpha = alpha, beta = beta, mar1 = mar1,
      mar2 = mar2, lower.tail = lower.tail)) 
}

"abvlog"<- 
function(x = 0.5, dep, plot = FALSE, add = FALSE, lty = 1, lwd = 1, col = 1,
         blty = 3, xlim = c(0,1), ylim = c(0.5,1), xlab = "", ylab = "", ...)
{
    if(mode(x) != "numeric" || any(x < 0,na.rm=TRUE) || any(x > 1,na.rm=TRUE))
        stop("invalid argument for `x'")
    if(length(dep) != 1 || mode(dep) != "numeric" || dep <= 0 ||
        dep > 1) stop("invalid argument for `dep'")
    if(plot || add) x <- seq(0, 1, length = 100)
    idep <- 1/dep
    a <- (x^idep + (1-x)^idep)^dep
    if(plot || add) {
        bvdepfn(x, a, add, lty, lwd, col, blty, xlab, ylab, xlim, ylim, ...)  
        return(invisible(list(x = x, y = a)))
    }
    a
}

"abvalog"<- 
function(x = 0.5, dep, asy = c(1,1), plot = FALSE, add = FALSE,
         lty = 1, lwd = 1, col = 1, blty = 3, xlim = c(0,1),
         ylim = c(0.5,1), xlab = "", ylab = "", ...)
{
    if(mode(x) != "numeric" || any(x < 0,na.rm=TRUE) || any(x > 1,na.rm=TRUE))
        stop("invalid argument for `x'")
    if(length(dep) != 1 || mode(dep) != "numeric" || dep <= 0 ||
        dep > 1) stop("invalid argument for `dep'")
    if(length(asy) != 2 || mode(asy) != "numeric" || min(asy) < 0 ||
       max(asy) > 1) stop("invalid argument for `asy'")
    if(plot || add) x <- seq(0, 1, length = 100)
    idep <- 1/dep
    a <- ((asy[1]*x)^idep + (asy[2]*(1-x))^idep)^dep +
        (1-asy[1])*x + (1-asy[2])*(1-x)    
    if(plot || add) {
        bvdepfn(x, a, add, lty, lwd, col, blty, xlab, ylab, xlim, ylim, ...) 
        return(invisible(list(x = x, y = a)))
    }
    a
}

"abvhr" <-
function(x = 0.5, dep, plot = FALSE, add = FALSE, lty = 1, lwd = 1, col = 1,
         blty = 3, xlim = c(0,1), ylim = c(0.5,1), xlab = "", ylab = "", ...)
{
    if(mode(x) != "numeric" || any(x < 0,na.rm=TRUE) || any(x > 1,na.rm=TRUE))
        stop("invalid argument for `x'")
    if(length(dep) != 1 || mode(dep) != "numeric" || dep <= 0)
        stop("invalid argument for `dep'")
    if(plot || add) x <- seq(0, 1, length = 100)
    fn <- function(z) z*pnorm(1/dep + dep * log(z/(1-z)) / 2)
    a <- fn(x) + fn(1-x)
    if(plot || add) {
        bvdepfn(x, a, add, lty, lwd, col, blty, xlab, ylab, xlim, ylim, ...)  
        return(invisible(list(x = x, y = a)))
    }
    a
}

"abvneglog"<- 
function(x = 0.5, dep, plot = FALSE, add = FALSE, lty = 1, lwd = 1, col = 1,
         blty = 3, xlim = c(0,1), ylim = c(0.5,1), xlab = "", ylab = "", ...)
{
    if(mode(x) != "numeric" || any(x < 0,na.rm=TRUE) || any(x > 1,na.rm=TRUE))
        stop("invalid argument for `x'")  
    if(length(dep) != 1 || mode(dep) != "numeric" || dep <= 0)
        stop("invalid argument for `dep'")
    if(plot || add) x <- seq(0, 1, length = 100)
    a <- 1 - (x^(-dep) + (1-x)^(-dep))^(-1/dep)
    if(plot || add) {
        bvdepfn(x, a, add, lty, lwd, col, blty, xlab, ylab, xlim, ylim, ...) 
        return(invisible(list(x = x, y = a)))
    }
    a
}

"abvaneglog"<- 
function(x = 0.5, dep, asy = c(1,1), plot = FALSE, add = FALSE,
         lty = 1, lwd = 1, col = 1, blty = 3, xlim = c(0,1),
         ylim = c(0.5,1), xlab = "", ylab = "", ...)
{
    if(mode(x) != "numeric" || any(x < 0,na.rm=TRUE) || any(x > 1,na.rm=TRUE))
        stop("invalid argument for `x'")
    if(length(dep) != 1 || mode(dep) != "numeric" || dep <= 0)
        stop("invalid argument for `dep'")
    if(plot || add) x <- seq(0, 1, length = 100)
    if(length(asy) != 2 || mode(asy) != "numeric" || min(asy) < 0 ||
       max(asy) > 1) stop("invalid argument for `asy'")
    a <- 1 - ((asy[1]*x)^(-dep) + (asy[2]*(1-x))^(-dep))^(-1/dep)
    if(plot || add) {
        bvdepfn(x, a, add, lty, lwd, col, blty, xlab, ylab, xlim, ylim, ...)  
        return(invisible(list(x = x, y = a)))
    }
    a
}

"abvbilog"<- 
function(x = 0.5, alpha, beta, plot = FALSE, add = FALSE,
         lty = 1, lwd = 1, col = 1, blty = 3, xlim = c(0,1),
         ylim = c(0.5,1), xlab = "", ylab = "", ...)
{
    if(mode(x) != "numeric" || any(x < 0,na.rm=TRUE) || any(x > 1,na.rm=TRUE))
        stop("invalid argument for `x'")
    if(length(alpha) != 1 || mode(alpha) != "numeric")
        stop("invalid argument for `alpha'")
    if(length(beta) != 1 || mode(beta) != "numeric")
        stop("invalid argument for `beta'")
    if(any(c(alpha,beta) <= 0) || any(c(alpha,beta) >= 1))
        stop("`alpha' and `beta' must be in the open interval (0,1)")
    if(plot || add) x <- seq(0, 1, length = 100)
    gma <- numeric(length(x))
    for(i in 1:length(x)) {
        gmafn <- function(z)
            (1-alpha) * x[i] * (1-z)^beta - (1-beta) * (1-x[i]) * z^alpha
        if(is.na(x[i])) gma[i] <- NA
        else if(x[i] == 0) gma[i] <- 0
        else if(x[i] == 1) gma[i] <- 1
        else gma[i] <- uniroot(gmafn, lower = 0, upper = 1,
                               tol = .Machine$double.eps^0.5)$root
    }
    a <- x * gma^(1-alpha) + (1-x) * (1 - gma)^(1-beta)
    if(plot || add) {
        bvdepfn(x, a, add, lty, lwd, col, blty, xlab, ylab, xlim, ylim, ...)  
        return(invisible(list(x = x, y = a)))
    }
    a
}

"abvnegbilog"<- 
function(x = 0.5, alpha, beta, plot = FALSE, add = FALSE,
         lty = 1, lwd = 1, col = 1, blty = 3, xlim = c(0,1),
         ylim = c(0.5,1), xlab = "", ylab = "", ...)
{
    if(mode(x) != "numeric" || any(x < 0,na.rm=TRUE) || any(x > 1,na.rm=TRUE))
        stop("invalid argument for `x'")
    if(length(alpha) != 1 || mode(alpha) != "numeric")
        stop("invalid argument for `alpha'")
    if(length(beta) != 1 || mode(beta) != "numeric")
        stop("invalid argument for `beta'")
    if(any(c(alpha,beta) <= 0))
        stop("`alpha' and `beta' must be non-negative")
    if(plot || add) x <- seq(0, 1, length = 100)
    gma <- numeric(length(x))
    for(i in 1:length(x)) {
        gmafn <- function(z)
            (1+alpha) * x[i] * z^alpha - (1+beta) * (1-x[i]) * (1-z)^beta
        if(is.na(x[i])) gma[i] <- NA
        else if(x[i] == 0) gma[i] <- 1
        else if(x[i] == 1) gma[i] <- 0
        else gma[i] <- uniroot(gmafn, lower = 0, upper = 1,
                               tol = .Machine$double.eps^0.5)$root
    }
    a <- 1 - x * gma^(1+alpha) - (1-x) * (1 - gma)^(1+beta)
    if(plot || add) {
        bvdepfn(x, a, add, lty, lwd, col, blty, xlab, ylab, xlim, ylim, ...)  
        return(invisible(list(x = x, y = a)))
    }
    a
}

"abvct" <-
function(x = 0.5, alpha, beta, plot = FALSE, add = FALSE,
         lty = 1, lwd = 1, col = 1, blty = 3, xlim = c(0,1),
         ylim = c(0.5,1), xlab = "", ylab = "", ...)
{
    if(mode(x) != "numeric" || any(x < 0,na.rm=TRUE) || any(x > 1,na.rm=TRUE))
        stop("invalid argument for `x'")
    if(length(alpha) != 1 || mode(alpha) != "numeric")
        stop("invalid argument for `alpha'")
    if(length(beta) != 1 || mode(beta) != "numeric")
        stop("invalid argument for `beta'")
    if(any(c(alpha,beta) <= 0))
        stop("`alpha' and `beta' must be non-negative")
    if(plot || add) x <- seq(0, 1, length = 100)
    u <- (alpha * (1-x)) / (alpha * (1-x) + beta * x)
    a <- (1-x) * pbeta(u, shape1 = alpha, shape2 = beta + 1) +
      x * pbeta(u, shape1 = alpha + 1, shape2 = beta, lower.tail = FALSE)
    if(plot || add) {
        bvdepfn(x, a, add, lty, lwd, col, blty, xlab, ylab, xlim, ylim, ...)  
        return(invisible(list(x = x, y = a)))
    }
    a
}

"abvpar" <-
function(x = 0.5, dep, asy = c(1,1), alpha, beta, model = c("log", "alog",
    "hr", "neglog", "aneglog", "bilog", "negbilog", "ct"),
     plot = FALSE, add = FALSE, lty = 1, lwd = 1, col = 1, blty = 3,
     xlim = c(0,1), ylim = c(0.5,1), xlab = "", ylab = "", ...)
{
  model <- match.arg(model)
  m1 <- c("bilog", "negbilog", "ct")
  m2 <- c(m1, "log", "hr", "neglog")
  m3 <- c("log", "alog", "hr", "neglog", "aneglog")
  if((model %in% m1) && !missing(dep))
    warning("ignoring `dep' argument")
  if((model %in% m2) && !missing(asy))
    warning("ignoring `asy' argument")
  if((model %in% m3) && !missing(alpha))
    warning("ignoring `alpha' argument")
  if((model %in% m3) && !missing(beta))
    warning("ignoring `beta' argument")
    
  switch(model,
    log = abvlog(x = x, dep = dep, plot = plot, add = add, lty = lty,
      lwd = lwd, col = col, blty = blty, xlim = xlim, ylim = ylim,
      xlab = xlab, ylab = ylab, ...),
    alog = abvalog(x = x, dep = dep, asy = asy, plot = plot, add = add,
      lty = lty, lwd = lwd, col = col, blty = blty, xlim = xlim, ylim = ylim,
      xlab = xlab, ylab = ylab, ...),
    hr = abvhr(x = x, dep = dep, plot = plot, add = add, lty = lty,
      lwd = lwd, col = col, blty = blty, xlim = xlim, ylim = ylim,
      xlab = xlab, ylab = ylab, ...),
    neglog = abvneglog(x = x, dep = dep, plot = plot, add = add, lty = lty,
      lwd = lwd, col = col, blty = blty, xlim = xlim, ylim = ylim,
      xlab = xlab, ylab = ylab, ...),
    aneglog = abvaneglog(x = x, dep = dep, asy = asy, plot = plot, add = add,
      lty = lty, lwd = lwd, col = col, blty = blty, xlim = xlim, ylim = ylim,
      xlab = xlab, ylab = ylab, ...),
    bilog = abvbilog(x = x, alpha = alpha, beta = beta, plot = plot,
      add = add, lty = lty, lwd = lwd, col = col, blty = blty, xlim = xlim,
      ylim = ylim, xlab = xlab, ylab = ylab, ...),
    negbilog = abvnegbilog(x = x, alpha = alpha, beta = beta, plot = plot,
      add = add, lty = lty, lwd = lwd, col = col, blty = blty, xlim = xlim,
      ylim = ylim, xlab = xlab, ylab = ylab, ...),
    ct = abvct(x = x, alpha = alpha, beta = beta, plot = plot, add = add,
      lty = lty, lwd = lwd, col = col, blty = blty, xlim = xlim, ylim = ylim,
      xlab = xlab, ylab = ylab, ...)) 
}

"abvnonpar"<- 
function(x = 0.5, data, nsloc1 = NULL, nsloc2 = NULL,
         method = c("cfg","pickands","deheuvels","hall","tdo"), convex = FALSE,
         wf = function(t) t, kmar = NULL, plot = FALSE, add = FALSE,
         lty = 1, lwd = 1, col = 1, blty = 3, xlim = c(0,1), ylim = c(0.5,1),
         xlab = "", ylab = "", ...)
{
    if(mode(x) != "numeric" || any(x < 0, na.rm=TRUE) ||
       any(x > 1, na.rm=TRUE)) stop("invalid argument for `x'")
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
        data <- mtransform(data, list(mle.m1, mle.m2))       
        # End transform
    }
    else {
        if(length(kmar) != 3 || mode(kmar) != "numeric")
            stop("`kmar' should be a numeric vector of length three")
        if(!is.null(nsloc1) || !is.null(nsloc2))
            warning("ignoring `nsloc1' and `nsloc2' arguments")
        data <- mtransform(data, kmar)
    }
    data <- na.omit(data)
    if(plot || add) x <- seq(0, 1, length = 100)
    method <- match.arg(method)
    mpmin <- function(a,b) {
      a[a > b] <- b[a > b]
      a
    }
    nn <- nrow(data)

    if(method == "pickands") {
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
    if(method == "deheuvels") {
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
        zvec <- sort(data[,1] / (data[,1] + data[,2]))
        if(any(table(zvec) > 1)) zvec <- jitter(zvec)
        zratio <- log(zvec) - log(1-zvec)
        qvec <- exp(cumsum(zratio) / nn)
        if(!convex) {
            step.pos <- as.numeric(cut(x, breaks = c(-0.1, zvec, 1))) - 1
            step.pos.r <- step.pos / nn
            a <- x^step.pos.r * (1-x)^(1 - step.pos.r) * qvec[nn]^wf(x) /
              c(rep(1, sum(step.pos == 0)), qvec[step.pos])
            a <- pmin(1, pmax(a, x, 1-x))
        }
        else {
            x2 <- seq(0, 1, length = 250)
            step.pos <- as.numeric(cut(x2, breaks = c(-0.1, zvec, 1))) - 1
            step.pos.r <- step.pos / nn
            a <- x2^step.pos.r * (1-x2)^(1 - step.pos.r) * qvec[nn]^wf(x2) /
              c(rep(1, sum(step.pos == 0)), qvec[step.pos])
            a <- pmin(1, pmax(a, x2, 1-x2))
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
    if(method == "hall") {
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
        bvdepfn(x, a, add, lty, lwd, col, blty, xlab, ylab, xlim, ylim, ...)  
        return(invisible(list(x = x, y = a)))
    }
    a
}

"dbvlog"<- 
function(x, dep, mar1 = c(0,1,0), mar2 = mar1, log = FALSE)
{
    if(length(dep) != 1 || mode(dep) != "numeric" || dep <= 0 ||
        dep > 1) stop("invalid argument for `dep'")
    if(is.null(dim(x))) dim(x) <- c(1,2)
    mar1 <- matrix(t(mar1), nrow = nrow(x), ncol = 3, byrow = TRUE)
    mar2 <- matrix(t(mar2), nrow = nrow(x), ncol = 3, byrow = TRUE)
    d <- numeric(nrow(x))
    x <- mtransform(x, list(mar1, mar2))
    ext <- apply(x,1,function(z) any(z %in% c(0,Inf)))
    d[ext] <- -Inf    
    if(any(!ext)) {
        x <- x[!ext, ,drop=FALSE]
        mar1 <- mar1[!ext, ,drop=FALSE]
        mar2 <- mar2[!ext, ,drop=FALSE]
        idep <- 1/dep
        z <- apply(x^idep,1,sum)^dep
        lx <- log(x)
        .expr1 <- (idep+mar1[,3])*lx[,1] + (idep+mar2[,3])*lx[,2] -
            log(mar1[,2]*mar2[,2])
        d[!ext] <- .expr1 + (1-2*idep)*log(z) + log(idep-1+z) - z
    }
    if(!log) d <- exp(d)
    d
}

"dbvalog"<- 
function(x, dep, asy = c(1,1), mar1 = c(0,1,0), mar2 = mar1, log = FALSE)
{
    if(length(dep) != 1 || mode(dep) != "numeric" || dep <= 0 ||
        dep > 1) stop("invalid argument for `dep'")
    if(length(asy) != 2 || mode(asy) != "numeric" || min(asy) < 0 ||
       max(asy) > 1) stop("invalid argument for `asy'")
    if(is.null(dim(x))) dim(x) <- c(1,2)
    mar1 <- matrix(t(mar1), nrow = nrow(x), ncol = 3, byrow = TRUE)
    mar2 <- matrix(t(mar2), nrow = nrow(x), ncol = 3, byrow = TRUE)
    d <- numeric(nrow(x))
    x <- mtransform(x, list(mar1, mar2))
    ext <- apply(x,1,function(z) any(z %in% c(0,Inf)))
    d[ext] <- -Inf
    if(any(!ext)) {
        x <- x[!ext, ,drop=FALSE]
        mar1 <- mar1[!ext, ,drop=FALSE]
        mar2 <- mar2[!ext, ,drop=FALSE]
        asy <- matrix(asy, ncol = 2, nrow = nrow(x), byrow = TRUE)
        idep <- 1/dep
        z <- apply((asy*x)^idep,1,sum)^dep
        v <- z + apply((1-asy)*x,1,sum)
        f1asy <- (idep)*log(asy)
        f2asy <- log(1-asy)
        lx <- log(x)
        fx <- (idep-1)*lx
        jac <- (1+mar1[,3])*lx[,1] + (1+mar2[,3])*lx[,2] -
            log(mar1[,2]*mar2[,2])
        .expr1 <- apply(f2asy,1,sum)
        .expr2 <- f2asy[,1] + f1asy[,2] + fx[,2]
        .expr3 <- f2asy[,2] + f1asy[,1] + fx[,1]
        .expr4 <- (1-idep)*log(z) + log(exp(.expr2)+exp(.expr3))
        .expr5 <- apply(cbind(f1asy,fx),1,sum) + (1-2*idep)*log(z) +
            log(idep-1+z)
        d[!ext] <- log(exp(.expr1)+exp(.expr4)+exp(.expr5))-v+jac
    }
    if(!log) d <- exp(d)
    d
}

"dbvhr" <-
function(x, dep, mar1 = c(0,1,0), mar2 = mar1, log = FALSE)
{
    if(length(dep) != 1 || mode(dep) != "numeric" || dep <= 0)
        stop("invalid argument for `dep'")
    if(is.null(dim(x))) dim(x) <- c(1,2)
    mar1 <- matrix(t(mar1), nrow = nrow(x), ncol = 3, byrow = TRUE)
    mar2 <- matrix(t(mar2), nrow = nrow(x), ncol = 3, byrow = TRUE)
    d <- numeric(nrow(x))
    x <- mtransform(x, list(mar1, mar2))
    ext <- apply(x,1,function(z) any(z %in% c(0,Inf)))
    d[ext] <- -Inf
    if(any(!ext)) {
        x <- x[!ext, ,drop=FALSE]
        mar1 <- mar1[!ext, ,drop=FALSE]
        mar2 <- mar2[!ext, ,drop=FALSE]
        fn <- function(x1, x2, nm = pnorm) x1 *
            nm(1/dep + dep * log(x1/x2) / 2)
        v <- fn(x[,1], x[,2]) + fn(x[,2], x[,1])
        .expr1 <- fn(x[,1], x[,2]) * fn(x[,2], x[,1]) +
            dep * fn(x[,1], x[,2], nm = dnorm) / 2
        lx <- log(x)
        jac <- mar1[,3]*lx[,1] + mar2[,3]*lx[,2] - log(mar1[,2]*mar2[,2])
        d[!ext] <- log(.expr1)+jac-v
    }
    if(!log) d <- exp(d)
    d    
}

"dbvneglog"<- 
function(x, dep, mar1 = c(0,1,0), mar2 = mar1, log = FALSE)
{
    if(length(dep) != 1 || mode(dep) != "numeric" || dep <= 0)
        stop("invalid argument for `dep'")
    if(is.null(dim(x))) dim(x) <- c(1,2)
    mar1 <- matrix(t(mar1), nrow = nrow(x), ncol = 3, byrow = TRUE)
    mar2 <- matrix(t(mar2), nrow = nrow(x), ncol = 3, byrow = TRUE)
    d <- numeric(nrow(x))
    x <- mtransform(x, list(mar1, mar2))
    ext <- apply(x,1,function(z) any(z %in% c(0,Inf)))
    d[ext] <- -Inf
    if(any(!ext)) {
        x <- x[!ext, ,drop=FALSE]
        mar1 <- mar1[!ext, ,drop=FALSE]
        mar2 <- mar2[!ext, ,drop=FALSE]
        idep <- 1/dep
        z <- apply(x^(-dep),1,sum)^(-idep)
        v <- apply(x,1,sum) - z
        lx <- log(x)
        fx <- (-dep-1)*lx
        jac <- (1+mar1[,3])*lx[,1] + (1+mar2[,3])*lx[,2] -
            log(mar1[,2]*mar2[,2])
        .expr1 <- (1+dep)*log(z) + log(exp(fx[,1])+exp(fx[,2]))
        .expr2 <- fx[,1] + fx[,2] + (1+2*dep)*log(z) + log(1+dep+z)
        d[!ext] <- log(1-exp(.expr1)+exp(.expr2))-v+jac
    }
    if(!log) d <- exp(d)
    d    
}

"dbvaneglog"<- 
function(x, dep, asy = c(1,1), mar1 = c(0,1,0), mar2 = mar1, log = FALSE)
{
    if(length(dep) != 1 || mode(dep) != "numeric" || dep <= 0)
        stop("invalid argument for `dep'")
    if(length(asy) != 2 || mode(asy) != "numeric" || min(asy) < 0 ||
       max(asy) > 1) stop("invalid argument for `asy'")
    if(is.null(dim(x))) dim(x) <- c(1,2)
    mar1 <- matrix(t(mar1), nrow = nrow(x), ncol = 3, byrow = TRUE)
    mar2 <- matrix(t(mar2), nrow = nrow(x), ncol = 3, byrow = TRUE)
    d <- numeric(nrow(x))
    x <- mtransform(x, list(mar1, mar2))
    ext <- apply(x,1,function(z) any(z %in% c(0,Inf)))
    d[ext] <- -Inf
    if(any(!ext)) {
        x <- x[!ext, ,drop=FALSE]
        mar1 <- mar1[!ext, ,drop=FALSE]
        mar2 <- mar2[!ext, ,drop=FALSE]
        asy <- matrix(asy, ncol = 2, nrow = nrow(x), byrow = TRUE)
        idep <- 1/dep
        z <- apply((asy*x)^(-dep),1,sum)^(-idep)
        v <- apply(x,1,sum) - z
        fasy <- (-dep)*log(asy)
        lx <- log(x)
        fx <- (-dep-1)*lx
        jac <- (1+mar1[,3])*lx[,1] + (1+mar2[,3])*lx[,2] -
            log(mar1[,2]*mar2[,2])
        .expr1 <- fasy[,1] + fx[,1]
        .expr2 <- fasy[,2] + fx[,2]
        .expr3 <- (1+dep)*log(z) + log(exp(.expr1)+exp(.expr2))
        .expr4 <- apply(cbind(fasy,fx),1,sum) + (1+2*dep)*log(z) + log(1+dep+z)
        d[!ext] <- log(1-exp(.expr3)+exp(.expr4))-v+jac
    }
    if(!log) d <- exp(d)
    d
}

"dbvbilog"<- 
function(x, alpha, beta, mar1 = c(0,1,0), mar2 = mar1, log = FALSE)
{
    if(length(alpha) != 1 || mode(alpha) != "numeric")
        stop("invalid argument for `alpha'")
    if(length(beta) != 1 || mode(beta) != "numeric")
        stop("invalid argument for `beta'")
    if(any(c(alpha,beta) <= 0) || any(c(alpha,beta) >= 1))
        stop("`alpha' and `beta' must be in the open interval (0,1)")
    if(is.null(dim(x))) dim(x) <- c(1,2)
    mar1 <- matrix(t(mar1), nrow = nrow(x), ncol = 3, byrow = TRUE)
    mar2 <- matrix(t(mar2), nrow = nrow(x), ncol = 3, byrow = TRUE)
    d <- numeric(nrow(x))
    x <- mtransform(x, list(mar1, mar2))
    ext <- apply(x,1,function(z) any(z %in% c(0,Inf)))
    d[ext] <- -Inf
    if(any(!ext)) {
        x <- x[!ext, ,drop=FALSE]
        mar1 <- mar1[!ext, ,drop=FALSE]
        mar2 <- mar2[!ext, ,drop=FALSE]
        gma <- numeric(nrow(x))
        for(i in 1:nrow(x)) {
            gmafn <- function(z)
                (1-alpha) * x[i,1] * (1-z)^beta - (1-beta) * x[i,2] * z^alpha
            if(any(is.na(x[i,]))) gma[i] <- NA
            else gma[i] <- uniroot(gmafn, lower = 0, upper = 1,
                               tol = .Machine$double.eps^0.5)$root
        }
        v <- x[,1] * gma^(1-alpha) + x[,2] * (1 - gma)^(1-beta)
        lx <- log(x)
        jac <- (1+mar1[,3])*lx[,1] + (1+mar2[,3])*lx[,2] -
            log(mar1[,2]*mar2[,2])
        .expr1 <- exp((1-alpha)*log(gma) + (1-beta)*log(1-gma))
        .expr2 <- exp(log(1-alpha) + log(beta) + (beta - 1)*log(1-gma) +
            lx[,1]) + exp(log(1-beta) + log(alpha) + (alpha - 1)*log(gma) +
            lx[,2])
        d[!ext] <- log(.expr1 + (1-alpha)*(1-beta)/.expr2) - v + jac
    }
    if(!log) d <- exp(d)
    d   
}

"dbvnegbilog"<- 
function(x, alpha, beta, mar1 = c(0,1,0), mar2 = mar1, log = FALSE)
{
    if(length(alpha) != 1 || mode(alpha) != "numeric")
        stop("invalid argument for `alpha'")
    if(length(beta) != 1 || mode(beta) != "numeric")
        stop("invalid argument for `beta'")
    if(any(c(alpha,beta) <= 0))
        stop("`alpha' and `beta' must be non-negative")
    if(is.null(dim(x))) dim(x) <- c(1,2)
    mar1 <- matrix(t(mar1), nrow = nrow(x), ncol = 3, byrow = TRUE)
    mar2 <- matrix(t(mar2), nrow = nrow(x), ncol = 3, byrow = TRUE)
    d <- numeric(nrow(x))
    x <- mtransform(x, list(mar1, mar2))
    ext <- apply(x,1,function(z) any(z %in% c(0,Inf)))
    d[ext] <- -Inf
    if(any(!ext)) {
        x <- x[!ext, ,drop=FALSE]
        mar1 <- mar1[!ext, ,drop=FALSE]
        mar2 <- mar2[!ext, ,drop=FALSE]
        gma <- numeric(nrow(x))
        for(i in 1:nrow(x)) {
            gmafn <- function(z)
                (1+alpha) * x[i,1] * z^alpha - (1+beta) * x[i,2] * (1-z)^beta
            if(any(is.na(x[i,]))) gma[i] <- NA
            else gma[i] <- uniroot(gmafn, lower = 0, upper = 1,
                               tol = .Machine$double.eps^0.5)$root
        }
        v <- x[,1] + x[,2] - x[,1] * gma^(1+alpha) - x[,2] * (1 - gma)^(1+beta)
        lx <- log(x)
        jac <- (1+mar1[,3])*lx[,1] + (1+mar2[,3])*lx[,2] -
            log(mar1[,2]*mar2[,2])
        .expr1 <- (1-gma^(1+alpha)) * (1 - (1-gma)^(1+beta))
        .expr2 <- exp(log(1+alpha) + log(1+beta) + alpha*log(gma) +
                      beta*log(1-gma))
        .expr3 <- exp(log(1+alpha) + log(alpha) + (alpha - 1)*log(gma) +
            lx[,1]) + exp(log(1+beta) + log(beta) + (beta - 1)*log(1-gma) +
            lx[,2])
        d[!ext] <- log(.expr1 + .expr2/.expr3) - v + jac
    }
    if(!log) d <- exp(d)
    d   
}

"dbvct"<- 
function(x, alpha, beta, mar1 = c(0,1,0), mar2 = mar1, log = FALSE)
{
    if(length(alpha) != 1 || mode(alpha) != "numeric")
        stop("invalid argument for `alpha'")
    if(length(beta) != 1 || mode(beta) != "numeric")
        stop("invalid argument for `beta'")
    if(any(c(alpha,beta) <= 0))
        stop("`alpha' and `beta' must be non-negative")
    if(is.null(dim(x))) dim(x) <- c(1,2)
    mar1 <- matrix(t(mar1), nrow = nrow(x), ncol = 3, byrow = TRUE)
    mar2 <- matrix(t(mar2), nrow = nrow(x), ncol = 3, byrow = TRUE)
    d <- numeric(nrow(x))
    x <- mtransform(x, list(mar1, mar2))
    ext <- apply(x,1,function(z) any(z %in% c(0,Inf)))
    d[ext] <- -Inf
    if(any(!ext)) {
        x <- x[!ext, ,drop=FALSE]
        mar1 <- mar1[!ext, ,drop=FALSE]
        mar2 <- mar2[!ext, ,drop=FALSE]
        u <- (alpha * x[,2]) / (alpha * x[,2] + beta * x[,1]) 
        v <- x[,2] * pbeta(u, shape1 = alpha, shape2 = beta + 1) + x[,1] *
          pbeta(u, shape1 = alpha + 1, shape2 = beta, lower.tail = FALSE)
        lx <- log(x)
        jac <- (1+mar1[,3])*lx[,1] + (1+mar2[,3])*lx[,2] -
            log(mar1[,2]*mar2[,2])
        .c1 <- alpha * beta / (alpha + beta + 1)
        .expr1 <- pbeta(u, shape1 = alpha, shape2 = beta + 1) *
          pbeta(u, shape1 = alpha + 1, shape2 = beta, lower.tail = FALSE)
        .expr2 <- dbeta(u, shape1 = alpha + 1, shape2 = beta + 1) /
          (alpha * x[,2] + beta * x[,1]) 
        d[!ext] <- log(.expr1 + .c1 * .expr2) - v + jac
    }
    if(!log) d <- exp(d)
    d   
}

"dbvevd" <-
function(x, dep, asy = c(1,1), alpha, beta, model = c("log", "alog",
    "hr", "neglog", "aneglog", "bilog", "negbilog", "ct"),
    mar1 = c(0,1,0), mar2 = mar1, log = FALSE)
{
  model <- match.arg(model)
  m1 <- c("bilog", "negbilog", "ct")
  m2 <- c(m1, "log", "hr", "neglog")
  m3 <- c("log", "alog", "hr", "neglog", "aneglog")
  if((model %in% m1) && !missing(dep))
    warning("ignoring `dep' argument")
  if((model %in% m2) && !missing(asy))
    warning("ignoring `asy' argument")
  if((model %in% m3) && !missing(alpha))
    warning("ignoring `alpha' argument")
  if((model %in% m3) && !missing(beta))
    warning("ignoring `beta' argument")
    
  switch(model,
    log = dbvlog(x = x, dep = dep, mar1 = mar1, mar2 = mar2, log = log),
    alog = dbvalog(x = x, dep = dep, asy = asy, mar1 = mar1,
      mar2 = mar2, log = log),
    hr = dbvhr(x = x, dep = dep, mar1 = mar1, mar2 = mar2, log = log),
    neglog = dbvneglog(x = x, dep = dep, mar1 = mar1, mar2 = mar2, log = log),
    aneglog = dbvaneglog(x = x, dep = dep, asy = asy, mar1 = mar1,
      mar2 = mar2, log = log),
    bilog = dbvbilog(x = x, alpha = alpha, beta = beta, mar1 = mar1,
      mar2 = mar2, log = log),
    negbilog = dbvnegbilog(x = x, alpha = alpha, beta = beta, mar1 = mar1,
      mar2 = mar2, log = log),
    ct = dbvct(x = x, alpha = alpha, beta = beta, mar1 = mar1,
      mar2 = mar2, log = log)) 
}

"fbvlog"<- 
function(x, start, ..., nsloc1 = NULL, nsloc2 = NULL, cshape = cscale, cscale = cloc, cloc = FALSE, std.err = TRUE, dsm = TRUE, corr = FALSE, method = "BFGS", warn.inf = TRUE)
{
    nlbvlog <- function(loc1, scale1, shape1, loc2, scale2, shape2, dep)
    {
        if(cshape) shape2 <- shape1
        if(cscale) scale2 <- scale1
        
        if(any(c(scale1,scale2) < 0.01) || dep < 0.1 || dep > 1)
            return(1e6)
        if(!is.null(nsloc1)) {
            ns <- numeric(length(loc.param1))
            for(i in 1:length(ns))
                ns[i] <- get(loc.param1[i])
            loc1 <- drop(nslocmat1 %*% ns)
        }
        else loc1 <- rep(loc1, length.out = nrow(x))
        if(cloc) loc2 <- loc1 else {
          if(!is.null(nsloc2)) {
            ns <- numeric(length(loc.param2))
            for(i in 1:length(ns))
                ns[i] <- get(loc.param2[i])
            loc2 <- drop(nslocmat2 %*% ns)
          }
          else loc2 <- rep(loc2, length.out = nrow(x))
        }
        if(spx$n.m1)
            m1l <- .C("nlgev", spx$x.m1, spx$n.m1, loc1[spx$na == 2], scale1,
                      shape1, dns = double(1), PACKAGE = "evd")$dns
        else m1l <- 0
        if(spx$n.m2)
            m2l <- .C("nlgev", spx$x.m2, spx$n.m2, loc2[spx$na == 1], scale2,
                      shape2, dns = double(1), PACKAGE = "evd")$dns
        else m2l <- 0   
        bvl <- .C("nlbvlog", spx$x1, spx$x2, spx$n, dep, loc1[spx$na == 0],
                  scale1, shape1, loc2[spx$na == 0], scale2, shape2,
                  dns = double(1), PACKAGE = "evd")$dns
        if(any(is.nan(c(m1l,m2l,bvl)))) {
            warning("NaN returned in likelihood")
            return(1e6)
        }
        if(any(c(m1l,m2l,bvl) == 1e6)) return(1e6)
        else return(m1l + m2l + bvl)
    }
    if(cloc && !identical(nsloc1, nsloc2))
      stop("nsloc1 and nsloc2 must be identical")
    if(!is.null(nsloc1)) {
        nsloc1 <- nsloc.transform(x, nsloc1)
        nslocmat1 <- cbind(1,as.matrix(nsloc1))
    }
    if(!is.null(nsloc2)) {
        nsloc2 <- nsloc.transform(x, nsloc2)
        nslocmat2 <- cbind(1,as.matrix(nsloc2))
    }
    loc.param1 <- paste("loc1", c("",names(nsloc1)), sep="")
    loc.param2 <- paste("loc2", c("",names(nsloc2)), sep="")
    param <- c(loc.param1, "scale1", "shape1")
    if(!cloc) param <- c(param, loc.param2) else loc.param2 <- NULL
    if(!cscale) param <- c(param, "scale2")
    if(!cshape) param <- c(param, "shape2")
    param <- c(param, "dep")
    nmdots <- names(list(...))
    start <- bvstart.vals(x, start, nsloc1, nsloc2, nmdots, param,
                          loc.param1, loc.param2, model = "log") 
    spx <- sep.bvdata(x)
    nm <- names(start)
    l <- length(nm)
    fixed.param <- list(...)[nmdots %in% param]
    if(any(!(param %in% c(nm,names(fixed.param)))))
        stop("unspecified parameters")
    prind <- (5:7)[c(!cscale, !cshape, TRUE)]
    f <- c(as.list(numeric(length(loc.param1))), formals(nlbvlog)[2:3],
           as.list(numeric(length(loc.param2))), formals(nlbvlog)[prind])
    names(f) <- param
    m <- match(nm, param)
    if(any(is.na(m))) 
        stop("`start' specifies unknown arguments")    
    formals(nlbvlog) <- c(f[m], f[-m])
    nllh <- function(p, ...) nlbvlog(p, ...)
    if(l > 1)
        body(nllh) <- parse(text = paste("nlbvlog(", paste("p[",1:l,
            "]", collapse = ", "), ", ...)"))
    start.arg <- c(list(p = unlist(start)), fixed.param)
    if(warn.inf && do.call("nllh", start.arg) == 1e6)
        warning("negative log-likelihood is infinite at starting values")
    opt <- optim(start, nllh, hessian = TRUE, ..., method = method)
    bvpost.optim(x, opt, nm, nsloc1, nsloc2, fixed.param, std.err, dsm, corr,
                 sym = FALSE, cmar = c(cloc, cscale, cshape), model = "log") 
}

"fbvalog"<- 
function(x, start, ..., sym = FALSE, nsloc1 = NULL, nsloc2 = NULL, cshape = cscale, cscale = cloc, cloc = FALSE, std.err = TRUE, dsm = TRUE, corr = FALSE, method = "BFGS", warn.inf = TRUE)
{   
    nlbvalog <- function(loc1, scale1, shape1, loc2, scale2, shape2,
                          asy1, asy2, dep)
    {
        if(sym) asy2 <- asy1
        if(cshape) shape2 <- shape1
        if(cscale) scale2 <- scale1
          
        if(any(c(scale1,scale2) < 0.01) || any(c(dep,asy1,asy2) > 1) ||
           any(c(asy1,asy2) < 0.001) || dep < 0.1)
            return(1e6)
        if(!is.null(nsloc1)) {
            ns <- numeric(length(loc.param1))
            for(i in 1:length(ns))
                ns[i] <- get(loc.param1[i])
            loc1 <- drop(nslocmat1 %*% ns)
        }
        else loc1 <- rep(loc1, length.out = nrow(x))
        if(cloc) loc2 <- loc1 else {
          if(!is.null(nsloc2)) {
            ns <- numeric(length(loc.param2))
            for(i in 1:length(ns))
              ns[i] <- get(loc.param2[i])
            loc2 <- drop(nslocmat2 %*% ns)
          }
          else loc2 <- rep(loc2, length.out = nrow(x))
        }
        if(spx$n.m1)
            m1l <- .C("nlgev", spx$x.m1, spx$n.m1, loc1[spx$na == 2], scale1,
                shape1, dns = double(1), PACKAGE = "evd")$dns
        else m1l <- 0
        if(spx$n.m2)
            m2l <- .C("nlgev", spx$x.m2, spx$n.m2, loc2[spx$na == 1], scale2,
                shape2, dns = double(1), PACKAGE = "evd")$dns
        else m2l <- 0
        bvl <- .C("nlbvalog", spx$x1, spx$x2, spx$n, dep, asy1, asy2,
                loc1[spx$na == 0], scale1, shape1, loc2[spx$na == 0],
                scale2, shape2, dns = double(1), PACKAGE = "evd")$dns
        if(any(is.nan(c(m1l,m2l,bvl)))) {
            warning("NaN returned in likelihood")
            return(1e6)
        }
        if(any(c(m1l,m2l,bvl) == 1e6)) return(1e6)
        else return(m1l + m2l + bvl)  
    }
    if(cloc && !identical(nsloc1, nsloc2))
      stop("nsloc1 and nsloc2 must be identical")
    if(!is.null(nsloc1)) {
        nsloc1 <- nsloc.transform(x, nsloc1)
        nslocmat1 <- cbind(1,as.matrix(nsloc1))
    }
    if(!is.null(nsloc2)) {
        nsloc2 <- nsloc.transform(x, nsloc2)
        nslocmat2 <- cbind(1,as.matrix(nsloc2))
    }
    loc.param1 <- paste("loc1", c("",names(nsloc1)), sep="")
    loc.param2 <- paste("loc2", c("",names(nsloc2)), sep="")
    param <- c(loc.param1, "scale1", "shape1")
    if(!cloc) param <- c(param, loc.param2) else loc.param2 <- NULL
    if(!cscale) param <- c(param, "scale2")
    if(!cshape) param <- c(param, "shape2")
    if(!sym) param <- c(param, "asy1", "asy2", "dep")
    else param <- c(param, "asy1", "dep")
    nmdots <- names(list(...))
    start <- bvstart.vals(x, start, nsloc1, nsloc2, nmdots, param, loc.param1,
      loc.param2, model = "alog")
    spx <- sep.bvdata(x)
    nm <- names(start)
    l <- length(nm)
    fixed.param <- list(...)[nmdots %in% param]
    if(any(!(param %in% c(nm,names(fixed.param)))))
        stop("unspecified parameters")
    prind <- (5:9)[c(!cscale, !cshape, TRUE, !sym, TRUE)]
    f <- c(as.list(numeric(length(loc.param1))), formals(nlbvalog)[2:3],
      as.list(numeric(length(loc.param2))), formals(nlbvalog)[prind])
    names(f) <- param
    m <- match(nm, param)
    if(any(is.na(m))) 
        stop("`start' specifies unknown arguments")
    formals(nlbvalog) <- c(f[m], f[-m])
    nllh <- function(p, ...) nlbvalog(p, ...)
    if(l > 1) body(nllh) <- parse(text = paste("nlbvalog(", paste("p[",1:l,"]",
      collapse = ", "), ", ...)"))
    start.arg <- c(list(p = unlist(start)), fixed.param)
    if(warn.inf && do.call("nllh", start.arg) == 1e6)
        warning("negative log-likelihood is infinite at starting values")
    opt <- optim(start, nllh, hessian = TRUE, ..., method = method)
    bvpost.optim(x, opt, nm, nsloc1, nsloc2, fixed.param, std.err, dsm, corr,
      sym = sym, cmar = c(cloc, cscale, cshape), model = "alog")
}

"fbvhr"<- 
function(x, start, ..., nsloc1 = NULL, nsloc2 = NULL, cshape = cscale, cscale = cloc, cloc = FALSE, std.err = TRUE, dsm = TRUE, corr = FALSE, method = "BFGS", warn.inf = TRUE)
{
    nlbvhr <- function(loc1, scale1, shape1, loc2, scale2, shape2,
                       dep)
    {
        if(cshape) shape2 <- shape1
        if(cscale) scale2 <- scale1
        
        if(any(c(scale1,scale2) < 0.01) || dep < 0.2 || dep > 10)
            return(1e6)
        if(!is.null(nsloc1)) {
            ns <- numeric(length(loc.param1))
            for(i in 1:length(ns))
                ns[i] <- get(loc.param1[i])
            loc1 <- drop(nslocmat1 %*% ns)
        }
        else loc1 <- rep(loc1, length.out = nrow(x))
        if(cloc) loc2 <- loc1 else {
          if(!is.null(nsloc2)) {
            ns <- numeric(length(loc.param2))
            for(i in 1:length(ns))
                ns[i] <- get(loc.param2[i])
            loc2 <- drop(nslocmat2 %*% ns)
          }
          else loc2 <- rep(loc2, length.out = nrow(x))
        }
        if(spx$n.m1)
            m1l <- .C("nlgev", spx$x.m1, spx$n.m1, loc1[spx$na == 2], scale1,
                shape1, dns = double(1), PACKAGE = "evd")$dns
        else m1l <- 0
        if(spx$n.m2)
            m2l <- .C("nlgev", spx$x.m2, spx$n.m2, loc2[spx$na == 1], scale2,
                shape2, dns = double(1), PACKAGE = "evd")$dns
        else m2l <- 0
        bvl <- .C("nlbvhr", spx$x1, spx$x2, spx$n, dep, loc1[spx$na == 0],
                scale1, shape1, loc2[spx$na == 0], scale2, shape2,
                dns = double(1), PACKAGE = "evd")$dns
        if(any(is.nan(c(m1l,m2l,bvl)))) {
            warning("NaN returned in likelihood")
            return(1e6)
        }
        if(any(c(m1l,m2l,bvl) == 1e6)) return(1e6)
        else return(m1l + m2l + bvl)
    }
    if(cloc && !identical(nsloc1, nsloc2))
      stop("nsloc1 and nsloc2 must be identical")
    if(!is.null(nsloc1)) {
        nsloc1 <- nsloc.transform(x, nsloc1)
        nslocmat1 <- cbind(1,as.matrix(nsloc1))
    }
    if(!is.null(nsloc2)) {
        nsloc2 <- nsloc.transform(x, nsloc2)
        nslocmat2 <- cbind(1,as.matrix(nsloc2))
    }
    loc.param1 <- paste("loc1", c("",names(nsloc1)), sep="")
    loc.param2 <- paste("loc2", c("",names(nsloc2)), sep="")
    param <- c(loc.param1, "scale1", "shape1")
    if(!cloc) param <- c(param, loc.param2) else loc.param2 <- NULL
    if(!cscale) param <- c(param, "scale2")
    if(!cshape) param <- c(param, "shape2")
    param <- c(param, "dep")
    nmdots <- names(list(...))
    start <- bvstart.vals(x, start, nsloc1, nsloc2, nmdots, param,
                          loc.param1, loc.param2, model = "hr") 
    spx <- sep.bvdata(x)
    nm <- names(start)
    l <- length(nm)
    fixed.param <- list(...)[nmdots %in% param]
    if(any(!(param %in% c(nm,names(fixed.param)))))
        stop("unspecified parameters")
    prind <- (5:7)[c(!cscale, !cshape, TRUE)]
    f <- c(as.list(numeric(length(loc.param1))), formals(nlbvhr)[2:3],
           as.list(numeric(length(loc.param2))), formals(nlbvhr)[prind])
    names(f) <- param
    m <- match(nm, param)
    if(any(is.na(m))) 
        stop("`start' specifies unknown arguments")    
    formals(nlbvhr) <- c(f[m], f[-m])
    nllh <- function(p, ...) nlbvhr(p, ...)
    if(l > 1)
        body(nllh) <- parse(text = paste("nlbvhr(", paste("p[",1:l,
            "]", collapse = ", "), ", ...)"))
    start.arg <- c(list(p = unlist(start)), fixed.param)
    if(warn.inf && do.call("nllh", start.arg) == 1e6)
        warning("negative log-likelihood is infinite at starting values")
    opt <- optim(start, nllh, hessian = TRUE, ..., method = method)
    bvpost.optim(x, opt, nm, nsloc1, nsloc2, fixed.param, std.err, dsm, corr,
      sym = FALSE, cmar = c(cloc, cscale, cshape), model = "hr")
}

"fbvneglog"<- 
function(x, start, ..., nsloc1 = NULL, nsloc2 = NULL, cshape = cscale, cscale = cloc, cloc = FALSE, std.err = TRUE, dsm = TRUE, corr = FALSE, method = "BFGS", warn.inf = TRUE)
{
    nlbvneglog <- function(loc1, scale1, shape1, loc2, scale2, shape2,
                           dep)
    {
        if(cshape) shape2 <- shape1
        if(cscale) scale2 <- scale1
        
        if(any(c(scale1,scale2) < 0.01) || dep < 0.05 || dep > 5)
            return(1e6)
        if(!is.null(nsloc1)) {
            ns <- numeric(length(loc.param1))
            for(i in 1:length(ns))
                ns[i] <- get(loc.param1[i])
            loc1 <- drop(nslocmat1 %*% ns)
        }
        else loc1 <- rep(loc1, length.out = nrow(x))
        if(cloc) loc2 <- loc1 else {
          if(!is.null(nsloc2)) {
            ns <- numeric(length(loc.param2))
            for(i in 1:length(ns))
                ns[i] <- get(loc.param2[i])
            loc2 <- drop(nslocmat2 %*% ns)
          }
          else loc2 <- rep(loc2, length.out = nrow(x))
        }
        if(spx$n.m1)
            m1l <- .C("nlgev", spx$x.m1, spx$n.m1, loc1[spx$na == 2], scale1,
                shape1, dns = double(1), PACKAGE = "evd")$dns
        else m1l <- 0
        if(spx$n.m2)
            m2l <- .C("nlgev", spx$x.m2, spx$n.m2, loc2[spx$na == 1], scale2,
                shape2, dns = double(1), PACKAGE = "evd")$dns
        else m2l <- 0
        bvl <- .C("nlbvneglog", spx$x1, spx$x2, spx$n, dep, loc1[spx$na == 0],
                  scale1, shape1, loc2[spx$na == 0], scale2, shape2,
                  dns = double(1), PACKAGE = "evd")$dns
        if(any(is.nan(c(m1l,m2l,bvl)))) {
            warning("NaN returned in likelihood")
            return(1e6)
        }
        if(any(c(m1l,m2l,bvl) == 1e6)) return(1e6)
        else return(m1l + m2l + bvl)
    }
    if(cloc && !identical(nsloc1, nsloc2))
      stop("nsloc1 and nsloc2 must be identical")
    if(!is.null(nsloc1)) {
        nsloc1 <- nsloc.transform(x, nsloc1)
        nslocmat1 <- cbind(1,as.matrix(nsloc1))
    }
    if(!is.null(nsloc2)) {
        nsloc2 <- nsloc.transform(x, nsloc2)
        nslocmat2 <- cbind(1,as.matrix(nsloc2))
    }
    loc.param1 <- paste("loc1", c("",names(nsloc1)), sep="")
    loc.param2 <- paste("loc2", c("",names(nsloc2)), sep="")
    param <- c(loc.param1, "scale1", "shape1")
    if(!cloc) param <- c(param, loc.param2) else loc.param2 <- NULL
    if(!cscale) param <- c(param, "scale2")
    if(!cshape) param <- c(param, "shape2")
    param <- c(param, "dep")
    nmdots <- names(list(...))
    start <- bvstart.vals(x, start, nsloc1, nsloc2, nmdots, param,
                          loc.param1, loc.param2, model = "neglog") 
    spx <- sep.bvdata(x)
    nm <- names(start)
    l <- length(nm)
    fixed.param <- list(...)[nmdots %in% param]
    if(any(!(param %in% c(nm,names(fixed.param)))))
        stop("unspecified parameters")
    prind <- (5:7)[c(!cscale, !cshape, TRUE)]
    f <- c(as.list(numeric(length(loc.param1))), formals(nlbvneglog)[2:3],
           as.list(numeric(length(loc.param2))), formals(nlbvneglog)[prind])
    names(f) <- param   
    m <- match(nm, param)
    if(any(is.na(m))) 
        stop("`start' specifies unknown arguments")    
    formals(nlbvneglog) <- c(f[m], f[-m])
    nllh <- function(p, ...) nlbvneglog(p, ...)
    if(l > 1)
        body(nllh) <- parse(text = paste("nlbvneglog(", paste("p[",1:l,
            "]", collapse = ", "), ", ...)"))
    start.arg <- c(list(p = unlist(start)), fixed.param)
    if(warn.inf && do.call("nllh", start.arg) == 1e6)
        warning("negative log-likelihood is infinite at starting values")
    opt <- optim(start, nllh, hessian = TRUE, ..., method = method)
    bvpost.optim(x, opt, nm, nsloc1, nsloc2, fixed.param, std.err, dsm, corr,
      sym = FALSE, cmar = c(cloc, cscale, cshape), model = "neglog")
}

"fbvaneglog"<- 
function(x, start, ..., sym = FALSE, nsloc1 = NULL, nsloc2 = NULL, cshape = cscale, cscale = cloc, cloc = FALSE, std.err = TRUE, dsm = TRUE, corr = FALSE, method = "BFGS", warn.inf = TRUE)
{
    nlbvaneglog <- function(loc1, scale1, shape1, loc2, scale2, shape2,
                            asy1, asy2, dep)
    {
        if(sym) asy2 <- asy1
        if(cshape) shape2 <- shape1
        if(cscale) scale2 <- scale1
        
        if(any(c(scale1,scale2) < 0.01) || any(c(asy1,asy2) > 1) ||
           any(c(asy1,asy2) < 0.001) || dep < 0.05 || dep > 5)
            return(1e6)
        if(!is.null(nsloc1)) {
            ns <- numeric(length(loc.param1))
            for(i in 1:length(ns))
                ns[i] <- get(loc.param1[i])
            loc1 <- drop(nslocmat1 %*% ns)
        }
        else loc1 <- rep(loc1, length.out = nrow(x))
        if(cloc) loc2 <- loc1 else {
          if(!is.null(nsloc2)) {
            ns <- numeric(length(loc.param2))
            for(i in 1:length(ns))
              ns[i] <- get(loc.param2[i])
            loc2 <- drop(nslocmat2 %*% ns)
          }
          else loc2 <- rep(loc2, length.out = nrow(x))
        }
        if(spx$n.m1)
            m1l <- .C("nlgev", spx$x.m1, spx$n.m1, loc1[spx$na == 2],
                scale1, shape1, dns = double(1), PACKAGE = "evd")$dns
        else m1l <- 0
        if(spx$n.m2)
            m2l <- .C("nlgev", spx$x.m2, spx$n.m2, loc2[spx$na == 1],
                scale2, shape2, dns = double(1), PACKAGE = "evd")$dns
        else m2l <- 0
        bvl <- .C("nlbvaneglog", spx$x1, spx$x2, spx$n, dep, asy1, asy2,
                loc1[spx$na == 0], scale1, shape1, loc2[spx$na == 0],
                scale2, shape2, dns = double(1), PACKAGE = "evd")$dns
        if(any(is.nan(c(m1l,m2l,bvl)))) {
            warning("NaN returned in likelihood")
            return(1e6)
        }
        if(any(c(m1l,m2l,bvl) == 1e6)) return(1e6)
        else return(m1l + m2l + bvl)
    }
    if(cloc && !identical(nsloc1, nsloc2))
      stop("nsloc1 and nsloc2 must be identical")
    if(!is.null(nsloc1)) {
        nsloc1 <- nsloc.transform(x, nsloc1)
        nslocmat1 <- cbind(1,as.matrix(nsloc1))
    }
    if(!is.null(nsloc2)) {
        nsloc2 <- nsloc.transform(x, nsloc2)
        nslocmat2 <- cbind(1,as.matrix(nsloc2))
    }
    loc.param1 <- paste("loc1", c("",names(nsloc1)), sep="")
    loc.param2 <- paste("loc2", c("",names(nsloc2)), sep="")
    param <- c(loc.param1, "scale1", "shape1")
    if(!cloc) param <- c(param, loc.param2) else loc.param2 <- NULL
    if(!cscale) param <- c(param, "scale2")
    if(!cshape) param <- c(param, "shape2")
    if(!sym) param <- c(param, "asy1", "asy2", "dep")
    else param <- c(param, "asy1", "dep")
    nmdots <- names(list(...))
    start <- bvstart.vals(x, start, nsloc1, nsloc2, nmdots, param,
      loc.param1, loc.param2, model = "aneglog") 
    spx <- sep.bvdata(x)
    nm <- names(start)
    l <- length(nm)
    fixed.param <- list(...)[nmdots %in% param]
    if(any(!(param %in% c(nm,names(fixed.param)))))
        stop("unspecified parameters")
    prind <- (5:9)[c(!cscale, !cshape, TRUE, !sym, TRUE)]
    f <- c(as.list(numeric(length(loc.param1))), formals(nlbvaneglog)[2:3],
      as.list(numeric(length(loc.param2))), formals(nlbvaneglog)[prind])
    names(f) <- param
    m <- match(nm, param)
    if(any(is.na(m))) 
        stop("`start' specifies unknown arguments")
    formals(nlbvaneglog) <- c(f[m], f[-m])
    nllh <- function(p, ...) nlbvaneglog(p, ...)
    if(l > 1) body(nllh) <- parse(text = paste("nlbvaneglog(",
      paste("p[",1:l,"]", collapse = ", "), ", ...)"))
    start.arg <- c(list(p = unlist(start)), fixed.param)
    if(warn.inf && do.call("nllh", start.arg) == 1e6)
        warning("negative log-likelihood is infinite at starting values")
    opt <- optim(start, nllh, hessian = TRUE, ..., method = method)
    bvpost.optim(x, opt, nm, nsloc1, nsloc2, fixed.param, std.err, dsm, corr,
      sym = sym, cmar = c(cloc, cscale, cshape), model = "aneglog")
}

"fbvbilog"<- 
function(x, start, ..., nsloc1 = NULL, nsloc2 = NULL, cshape = cscale, cscale = cloc, cloc = FALSE, std.err = TRUE, dsm = TRUE, corr = FALSE, method = "BFGS", warn.inf = TRUE)
{
    nlbvbilog <- function(loc1, scale1, shape1, loc2, scale2, shape2,
                          alpha, beta)
    {
        if(cshape) shape2 <- shape1
        if(cscale) scale2 <- scale1
        
        if(any(c(scale1,scale2) < 0.01) || any(c(alpha,beta) < 0.1) ||
           any(c(alpha,beta) > 0.999))
            return(1e6)
        if(!is.null(nsloc1)) {
            ns <- numeric(length(loc.param1))
            for(i in 1:length(ns))
                ns[i] <- get(loc.param1[i])
            loc1 <- drop(nslocmat1 %*% ns)
        }
        else loc1 <- rep(loc1, length.out = nrow(x))
        if(cloc) loc2 <- loc1 else {
          if(!is.null(nsloc2)) {
            ns <- numeric(length(loc.param2))
            for(i in 1:length(ns))
                ns[i] <- get(loc.param2[i])
            loc2 <- drop(nslocmat2 %*% ns)
          }
          else loc2 <- rep(loc2, length.out = nrow(x))
        }
        if(spx$n.m1)
            m1l <- .C("nlgev", spx$x.m1,  spx$n.m1, loc1[spx$na == 2],
                scale1, shape1, dns = double(1), PACKAGE = "evd")$dns
        else m1l <- 0
        if(spx$n.m2)
            m2l <- .C("nlgev", spx$x.m2,  spx$n.m2, loc2[spx$na == 1],
                scale2, shape2, dns = double(1), PACKAGE = "evd")$dns
        else m2l <- 0
        bvl <- .C("nlbvbilog", spx$x1, spx$x2, spx$n, alpha, beta,
                loc1[spx$na == 0], scale1, shape1, loc2[spx$na == 0],
                scale2, shape2, dns = double(1), PACKAGE = "evd")$dns
        if(any(is.nan(c(m1l,m2l,bvl)))) {
            warning("NaN returned in likelihood")
            return(1e6)
        }
        if(any(c(m1l,m2l,bvl) == 1e6)) return(1e6)
        else return(m1l + m2l + bvl)
    }
    if(cloc && !identical(nsloc1, nsloc2))
      stop("nsloc1 and nsloc2 must be identical")
    if(!is.null(nsloc1)) {
        nsloc1 <- nsloc.transform(x, nsloc1)
        nslocmat1 <- cbind(1,as.matrix(nsloc1))
    }
    if(!is.null(nsloc2)) {
        nsloc2 <- nsloc.transform(x, nsloc2)
        nslocmat2 <- cbind(1,as.matrix(nsloc2))
    }
    loc.param1 <- paste("loc1", c("",names(nsloc1)), sep="")
    loc.param2 <- paste("loc2", c("",names(nsloc2)), sep="")
    param <- c(loc.param1, "scale1", "shape1")
    if(!cloc) param <- c(param, loc.param2) else loc.param2 <- NULL
    if(!cscale) param <- c(param, "scale2")
    if(!cshape) param <- c(param, "shape2")
    param <- c(param, "alpha", "beta")
    nmdots <- names(list(...))
    start <- bvstart.vals(x, start, nsloc1, nsloc2, nmdots, param,
                          loc.param1, loc.param2, model = "bilog") 
    spx <- sep.bvdata(x)
    nm <- names(start)
    l <- length(nm)
    fixed.param <- list(...)[nmdots %in% param]
    if(any(!(param %in% c(nm,names(fixed.param)))))
        stop("unspecified parameters")
    prind <- (5:8)[c(!cscale, !cshape, TRUE, TRUE)]
    f <- c(as.list(numeric(length(loc.param1))), formals(nlbvbilog)[2:3],
           as.list(numeric(length(loc.param2))), formals(nlbvbilog)[prind])
    names(f) <- param
    m <- match(nm, param)
    if(any(is.na(m))) 
        stop("`start' specifies unknown arguments")    
    formals(nlbvbilog) <- c(f[m], f[-m])
    nllh <- function(p, ...) nlbvbilog(p, ...)
    if(l > 1)
        body(nllh) <- parse(text = paste("nlbvbilog(", paste("p[",1:l,
            "]", collapse = ", "), ", ...)"))
    start.arg <- c(list(p = unlist(start)), fixed.param)
    if(warn.inf && do.call("nllh", start.arg) == 1e6)
        warning("negative log-likelihood is infinite at starting values")
    opt <- optim(start, nllh, hessian = TRUE, ..., method = method)
    bvpost.optim(x, opt, nm, nsloc1, nsloc2, fixed.param, std.err, dsm, corr,
      sym = FALSE, cmar = c(cloc, cscale, cshape), model = "bilog")
}

"fbvnegbilog"<- 
function(x, start, ..., nsloc1 = NULL, nsloc2 = NULL, cshape = cscale, cscale = cloc, cloc = FALSE, std.err = TRUE, dsm = TRUE, corr = FALSE, method = "BFGS", warn.inf = TRUE)
{
    nlbvnegbilog <- function(loc1, scale1, shape1, loc2, scale2, shape2,
                             alpha, beta)
    {
        if(cshape) shape2 <- shape1
        if(cscale) scale2 <- scale1
        
        if(any(c(scale1,scale2) < 0.01) || any(c(alpha,beta) < 0.1) ||
           any(c(alpha,beta) > 20))
            return(1e6)
        if(!is.null(nsloc1)) {
            ns <- numeric(length(loc.param1))
            for(i in 1:length(ns))
                ns[i] <- get(loc.param1[i])
            loc1 <- drop(nslocmat1 %*% ns)
        }
        else loc1 <- rep(loc1, length.out = nrow(x))
        if(cloc) loc2 <- loc1 else {
          if(!is.null(nsloc2)) {
            ns <- numeric(length(loc.param2))
            for(i in 1:length(ns))
                ns[i] <- get(loc.param2[i])
            loc2 <- drop(nslocmat2 %*% ns)
          }
          else loc2 <- rep(loc2, length.out = nrow(x))
        }
        if(spx$n.m1)
            m1l <- .C("nlgev", spx$x.m1, spx$n.m1, loc1[spx$na == 2],
                scale1, shape1, dns = double(1), PACKAGE = "evd")$dns
        else m1l <- 0
        if(spx$n.m2)
            m2l <- .C("nlgev", spx$x.m2, spx$n.m2, loc2[spx$na == 1],
                scale2, shape2, dns = double(1), PACKAGE = "evd")$dns
        else m2l <- 0
        bvl <- .C("nlbvnegbilog", spx$x1, spx$x2, spx$n, alpha, beta,
                  loc1[spx$na == 0], scale1, shape1, loc2[spx$na == 0],
                  scale2, shape2, dns = double(1), PACKAGE = "evd")$dns
        if(any(is.nan(c(m1l,m2l,bvl)))) {
            warning("NaN returned in likelihood")
            return(1e6)
        }
        if(any(c(m1l,m2l,bvl) == 1e6)) return(1e6)
        else return(m1l + m2l + bvl)
    }
    if(cloc && !identical(nsloc1, nsloc2))
      stop("nsloc1 and nsloc2 must be identical")
    if(!is.null(nsloc1)) {
        nsloc1 <- nsloc.transform(x, nsloc1)
        nslocmat1 <- cbind(1,as.matrix(nsloc1))
    }
    if(!is.null(nsloc2)) {
        nsloc2 <- nsloc.transform(x, nsloc2)
        nslocmat2 <- cbind(1,as.matrix(nsloc2))
    }
    loc.param1 <- paste("loc1", c("",names(nsloc1)), sep="")
    loc.param2 <- paste("loc2", c("",names(nsloc2)), sep="")
    param <- c(loc.param1, "scale1", "shape1")
    if(!cloc) param <- c(param, loc.param2) else loc.param2 <- NULL
    if(!cscale) param <- c(param, "scale2")
    if(!cshape) param <- c(param, "shape2")
    param <- c(param, "alpha", "beta")
    nmdots <- names(list(...))
    start <- bvstart.vals(x, start, nsloc1, nsloc2, nmdots, param,
                          loc.param1, loc.param2, model = "negbilog") 
    spx <- sep.bvdata(x)
    nm <- names(start)
    l <- length(nm)
    fixed.param <- list(...)[nmdots %in% param]
    if(any(!(param %in% c(nm,names(fixed.param)))))
        stop("unspecified parameters")
    prind <- (5:8)[c(!cscale, !cshape, TRUE, TRUE)]
    f <- c(as.list(numeric(length(loc.param1))), formals(nlbvnegbilog)[2:3],
           as.list(numeric(length(loc.param2))), formals(nlbvnegbilog)[prind])
    names(f) <- param
    m <- match(nm, param)
    if(any(is.na(m))) 
        stop("`start' specifies unknown arguments")    
    formals(nlbvnegbilog) <- c(f[m], f[-m])
    nllh <- function(p, ...) nlbvnegbilog(p, ...)
    if(l > 1)
        body(nllh) <- parse(text = paste("nlbvnegbilog(", paste("p[",1:l,
            "]", collapse = ", "), ", ...)"))
    start.arg <- c(list(p = unlist(start)), fixed.param)
    if(warn.inf && do.call("nllh", start.arg) == 1e6)
        warning("negative log-likelihood is infinite at starting values") 
    opt <- optim(start, nllh, hessian = TRUE, ..., method = method)
    bvpost.optim(x, opt, nm, nsloc1, nsloc2, fixed.param, std.err, dsm, corr,
      sym = FALSE, cmar = c(cloc, cscale, cshape), model = "negbilog")
}

"fbvct"<- 
function(x, start, ..., sym = FALSE, nsloc1 = NULL, nsloc2 = NULL, cshape = cscale, cscale = cloc, cloc = FALSE, std.err = TRUE, dsm = TRUE, corr = FALSE, method = "BFGS", warn.inf = TRUE)
{
    nlbvct <- function(loc1, scale1, shape1, loc2, scale2, shape2,
                       alpha, beta)
    {
        if(sym) beta <- alpha
        if(cshape) shape2 <- shape1
        if(cscale) scale2 <- scale1
        
        if(any(c(scale1,scale2) < 0.01) || any(c(alpha,beta) < 0.001) ||
           any(c(alpha,beta) > 30))
            return(1e6)
        if(!is.null(nsloc1)) {
            ns <- numeric(length(loc.param1))
            for(i in 1:length(ns))
                ns[i] <- get(loc.param1[i])
            loc1 <- drop(nslocmat1 %*% ns)
        }
        else loc1 <- rep(loc1, length.out = nrow(x))
        if(cloc) loc2 <- loc1 else {
          if(!is.null(nsloc2)) {
            ns <- numeric(length(loc.param2))
            for(i in 1:length(ns))
              ns[i] <- get(loc.param2[i])
            loc2 <- drop(nslocmat2 %*% ns)
          }
          else loc2 <- rep(loc2, length.out = nrow(x))
        }
        if(spx$n.m1)
            m1l <- .C("nlgev", spx$x.m1, spx$n.m1, loc1[spx$na == 2],
                scale1, shape1, dns = double(1), PACKAGE = "evd")$dns
        else m1l <- 0
        if(spx$n.m2)
            m2l <- .C("nlgev", spx$x.m2, spx$n.m2, loc2[spx$na == 1],
                scale2, shape2, dns = double(1), PACKAGE = "evd")$dns
        else m2l <- 0
        bvl <- .C("nlbvct", spx$x1, spx$x2, spx$n, alpha, beta,
                loc1[spx$na == 0], scale1, shape1, loc2[spx$na == 0],
                scale2, shape2, dns = double(1), PACKAGE = "evd")$dns
        if(any(is.nan(c(m1l,m2l,bvl)))) {
            warning("NaN returned in likelihood")
            return(1e6)
        }
        if(any(c(m1l,m2l,bvl) == 1e6)) return(1e6)
        else return(m1l + m2l + bvl)
    }
    if(cloc && !identical(nsloc1, nsloc2))
      stop("nsloc1 and nsloc2 must be identical")
    if(!is.null(nsloc1)) {
        nsloc1 <- nsloc.transform(x, nsloc1)
        nslocmat1 <- cbind(1,as.matrix(nsloc1))
    }
    if(!is.null(nsloc2)) {
        nsloc2 <- nsloc.transform(x, nsloc2)
        nslocmat2 <- cbind(1,as.matrix(nsloc2))
    }
    loc.param1 <- paste("loc1", c("",names(nsloc1)), sep="")
    loc.param2 <- paste("loc2", c("",names(nsloc2)), sep="")
    param <- c(loc.param1, "scale1", "shape1")
    if(!cloc) param <- c(param, loc.param2) else loc.param2 <- NULL
    if(!cscale) param <- c(param, "scale2")
    if(!cshape) param <- c(param, "shape2")
    if(!sym) param <- c(param, "alpha", "beta")
    else param <- c(param, "alpha")
    nmdots <- names(list(...))
    start <- bvstart.vals(x, start, nsloc1, nsloc2, nmdots, param,
      loc.param1, loc.param2, model = "ct") 
    spx <- sep.bvdata(x)
    nm <- names(start)
    l <- length(nm)
    fixed.param <- list(...)[nmdots %in% param]
    if(any(!(param %in% c(nm,names(fixed.param)))))
        stop("unspecified parameters")
    prind <- (5:8)[c(!cscale, !cshape, TRUE, !sym)]
    f <- c(as.list(numeric(length(loc.param1))), formals(nlbvct)[2:3],
      as.list(numeric(length(loc.param2))), formals(nlbvct)[prind])
    names(f) <- param
    m <- match(nm, param)
    if(any(is.na(m))) 
        stop("`start' specifies unknown arguments")
    formals(nlbvct) <- c(f[m], f[-m])
    nllh <- function(p, ...) nlbvct(p, ...)
    if(l > 1) body(nllh) <- parse(text = paste("nlbvct(",
      paste("p[",1:l,"]", collapse = ", "), ", ...)"))
    start.arg <- c(list(p = unlist(start)), fixed.param)
    if(warn.inf && do.call("nllh", start.arg) == 1e6)
        warning("negative log-likelihood is infinite at starting values")
    opt <- optim(start, nllh, hessian = TRUE, ..., method = method)
    bvpost.optim(x, opt, nm, nsloc1, nsloc2, fixed.param, std.err, dsm, corr,
      sym = sym, cmar = c(cloc, cscale, cshape), model = "ct")
}

"fbvevd" <-
function(x, model = c("log", "alog", "hr", "neglog", "aneglog", "bilog", "negbilog", "ct"), start, ..., sym = FALSE, nsloc1 = NULL, nsloc2 = NULL, cshape = cscale, cscale = cloc, cloc = FALSE, std.err = TRUE, dsm = TRUE, corr = FALSE, method = "BFGS", warn.inf = TRUE)
{
  call <- match.call()
  model <- match.arg(model)
  if(sym && !(model %in% c("alog","aneglog","ct")))
    warning("Argument `sym' was ignored")
  
  ft <- switch(model,
    log = fbvlog(x=x, start=start, ..., nsloc1=nsloc1, nsloc2=nsloc2,
      cshape=cshape, cscale=cscale, cloc=cloc, std.err=std.err,
      dsm=dsm, corr=corr, method=method, warn.inf=warn.inf),
    alog = fbvalog(x=x, start=start, ..., sym=sym, nsloc1=nsloc1,
      nsloc2=nsloc2, cshape=cshape, cscale=cscale, cloc=cloc,
      std.err=std.err, dsm=dsm, corr=corr, method=method,
      warn.inf=warn.inf),
    hr = fbvhr(x=x, start=start, ..., nsloc1=nsloc1, nsloc2=nsloc2,
      cshape=cshape, cscale=cscale, cloc=cloc, std.err=std.err, dsm=dsm,
      corr=corr, method=method, warn.inf=warn.inf),
    neglog = fbvneglog(x=x, start=start, ..., nsloc1=nsloc1, nsloc2=nsloc2,
      cshape=cshape, cscale=cscale, cloc=cloc, std.err=std.err, dsm=dsm,
      corr=corr, method=method, warn.inf=warn.inf),
    aneglog = fbvaneglog(x=x, start=start, ..., sym=sym, nsloc1=nsloc1,
      nsloc2=nsloc2, cshape=cshape, cscale=cscale, cloc=cloc,
      std.err=std.err, dsm=dsm, corr=corr, method=method,
      warn.inf=warn.inf),
    bilog = fbvbilog(x=x, start=start, ..., nsloc1=nsloc1, nsloc2=nsloc2,
      cshape=cshape, cscale=cscale, cloc=cloc, std.err=std.err, dsm=dsm,
      corr=corr, method=method, warn.inf=warn.inf),
    negbilog = fbvnegbilog(x=x, start=start, ..., nsloc1=nsloc1,
      nsloc2=nsloc2, cshape=cshape, cscale=cscale, cloc=cloc,
      std.err=std.err, dsm=dsm, corr= corr, method=method,
      warn.inf=warn.inf),
    ct = fbvct(x=x, start=start, ..., sym=sym, nsloc1=nsloc1,
      nsloc2=nsloc2, cshape=cshape, cscale=cscale, cloc=cloc,
      std.err=std.err, dsm=dsm, corr=corr, method=method,
      warn.inf=warn.inf))

  structure(c(ft, call = call), class = c("bvevd","evd"))
}

"print.bvevd" <-  function(x, digits = max(3, getOption("digits") - 3), ...) 
{
    cat("\nCall:", deparse(x$call), "\n")
    cat("Deviance:", deviance(x), "\n")
    cat("AIC:", AIC(x), "\n")
    cat("\nEstimates\n")
    print.default(format(fitted(x), digits = digits), print.gap = 2, 
        quote = FALSE)
    if(!is.null(std.errors(x))) {
      cat("\nStandard Errors\n")
      print.default(format(std.errors(x), digits = digits),
          print.gap = 2, quote = FALSE)
    }
    if(!is.null(x$corr)) {
      cat("\nCorrelations\n")
      print.default(format(x$corr, digits = digits), print.gap = 2, 
          quote = FALSE)
    }
    if(!is.null(x$dep.summary)) {
      cat("\nDependence Structure\n")
      cat("  Dependence One:", x$dep.summary[1], "\n")
      cat("  Dependence Two:", x$dep.summary[2], "\n")
      cat("  Asymmetry:", x$dep.summary[3], "\n")
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

"plot.bvevd" <-  function(x, mar = 0, which = 1:4, main =
     c("Conditional Plot One", "Conditional Plot Two", "Density Plot",
       "Dependence Function"),
     ask = nb.fig < length(which) && dev.interactive(), ci = TRUE,
     jitter = FALSE, grid = 50, nplty = 2, blty = 3, method = "cfg",
     convex = FALSE, wf = function(t) t, ...) 
{
    if (!inherits(x, "bvevd")) 
        stop("Use only with `bvevd' objects")
    nb.fig <- prod(par("mfcol"))
    if(mar == 1 || mar == 2) {
        indx <- paste(c("loc","scale","shape"), as.character(mar), sep="")
        tdata <- na.omit(x$tdata[, mar])
        n <- length(tdata)
        param <- x$param[indx]
        names(param) <- c("loc","scale","shape")
        gev.mar <- structure(list(param = param, tdata = tdata, n = n,
           loc = param["loc"]), class = c("gev", "uvevd", "evd"))
        if(missing(which)) which <- 1:4
        if(missing(main)) main <- c("Probability Plot", "Quantile Plot",
           "Density Plot", "Return Level Plot")
        plot(gev.mar, which = which, main = main, ask = ask, ci = ci,
             jitter = jitter, ...)
        return(invisible(x))
    }
    if (!is.numeric(which) || any(which < 1) || any(which > 4)) 
        stop("`which' must be in 1:4")
    show <- rep(FALSE, 4)
    show[which] <- TRUE
    if (ask) {
        op <- par(ask = TRUE)
        on.exit(par(op))
    }
    if (show[1]) {
        bvcpp(x, mar = 1, ci = ci, main = main[1], xlim = c(0,1),
              ylim = c(0,1), ...)
    }
    if (show[2]) {
        bvcpp(x, mar = 2, ci = ci, main = main[2], xlim = c(0,1),
              ylim = c(0,1), ...)
    }
    if (show[3]) {
        bvdens(x, jitter = jitter, grid = grid, main = main[3], ...)
    }
    if (show[4]) {
        bvdp(x, nplty = nplty, blty = blty, method = method,
             wf = wf, main = main[4], ...)
    }
    invisible(x)
} 

"bvcpp" <-  function(x, mar = 2, ci = TRUE, main = "Conditional Probability Plot", xlab = "Empirical", ylab = "Model", ...)
{ 
    data <- x$tdata
    mle.m1 <- x$param[c("loc1","scale1","shape1")]
    mle.m2 <- x$param[c("loc2","scale2","shape2")]
    data <- exp(-mtransform(data, list(mle.m1, mle.m2)))
    data <- na.omit(data)
    n <- nrow(data)
    ppx <- ppoints(n)
    if(x$model %in% c("log","hr","neglog"))
        probs <- ccop(data[,1], data[,2], mar = mar, dep = x$param["dep"],
                      model = x$model)
    if(x$model  %in% c("alog","aneglog"))
        probs <- ccop(data[,1], data[,2], mar = mar, dep = x$param["dep"],
                      asy = x$param[c("asy1","asy2")], model = x$model)
    if(x$model  %in% c("bilog","negbilog","ct"))
        probs <- ccop(data[,1], data[,2], mar = mar, alpha = x$param["alpha"],
                      beta = x$param["beta"], model = x$model)
    probs <- sort(probs)
    if(!ci) {
        plot(ppx, probs, main = main, xlab = xlab, ylab = ylab, ...)
        abline(0, 1)
    }
    else {
        samp <- runif(n*99)
        samp <- matrix(samp, n, 99)
        samp <- apply(samp, 2, sort)
        samp <- apply(samp, 1, sort)
        env <- t(samp[c(3,97),])
        matplot(ppx, cbind(probs, env), main = main, xlab = xlab,
                ylab = ylab, type = "pnn", pch = 4, ...)
        xyuser <- par("usr")
        smidge <- min(diff(c(xyuser[1], ppx, xyuser[2])))/2
        segments(ppx-smidge, env[,1], ppx+smidge, env[,1])
        segments(ppx-smidge, env[,2], ppx+smidge, env[,2])
        abline(0, 1)
    }
    invisible(list(x = ppx, y = probs))
}

"ccop" <- function(x1, x2, mar, dep, asy, alpha, beta, model)
{
    model <- match(model, c("log","alog","hr","neglog","aneglog",
                            "bilog","negbilog","ct"))
    if(model <= 5) alpha <- beta <- 1
    else dep <- 1
    if(model != 2 && model != 5) asy <- c(1,1)
    n <- length(x1)
    .C("ccop",
            as.double(x1), as.double(x2), as.integer(mar), as.double(dep),
            as.double(asy[1]), as.double(asy[2]), as.double(alpha),
            as.double(beta), as.integer(n), as.integer(model),
            ccop = double(n), PACKAGE = "evd")$ccop
}

"bvdens" <-  function(x, jitter = FALSE, grid = 50, main = "Density Plot", xlab = "", ylab = "", ...)
{
    xlimit <- range(x$tdata[,1], na.rm = TRUE)
    ylimit <- range(x$tdata[,2], na.rm = TRUE)
    xlimit[1] <- xlimit[1] - diff(xlimit) * 0.1
    xlimit[2] <- xlimit[2] + diff(xlimit) * 0.1
    ylimit[1] <- ylimit[1] - diff(ylimit) * 0.1
    ylimit[2] <- ylimit[2] + diff(ylimit) * 0.1
    xvec <- seq(xlimit[1], xlimit[2], length = grid)
    yvec <- seq(ylimit[1], ylimit[2], length = grid)
    xyvals <- expand.grid(xvec, yvec)
    mar1 <- x$param[c("loc1","scale1","shape1")]
    mar2 <- x$param[c("loc2","scale2","shape2")]
    if(x$model %in% c("log","hr","neglog"))
        dfunargs <- list(dep = x$param["dep"], mar1 = mar1, mar2 = mar2)
    if(x$model  %in% c("alog","aneglog"))
        dfunargs <- list(dep = x$param["dep"],
            asy = x$param[c("asy1","asy2")], mar1 = mar1, mar2 = mar2)
    if(x$model  %in% c("bilog","negbilog","ct"))
        dfunargs <- list(alpha = x$param["alpha"], beta = x$param["beta"],
            mar1 = mar1, mar2 = mar2)
    dfunargs <- c(list(x = xyvals, model = x$model), dfunargs)
    dens <- do.call("dbvevd", dfunargs)
    dens <- matrix(dens, nrow = grid, ncol = grid)
    contour(xvec, yvec, dens, main = main, xlab = xlab, ylab = ylab, ...)
    data <- na.omit(x$tdata)
    if(jitter) {
        data[,1] <- jitter(data[,1])
        data[,2] <- jitter(data[,2])
    }
    points(data, pch = 4)
    invisible(list(x = xyvals, y = dens))
}

"bvdp" <- function(x, method = "cfg", convex = FALSE, wf = function(t) t, add = FALSE, lty = 1, nplty = 2, blty = 3, main = "Dependence Function", xlab = "", ylab = "", ...)
{
    abvnonpar(data = x$data, nsloc1 = x$nsloc1, nsloc2 = x$nsloc2,
              method = method, convex = convex, wf = wf,
              plot = TRUE, lty = nplty, blty = blty, main = main,
              xlab = xlab, ylab = ylab, add = add, ...)
    if(x$model %in% c("log","hr","neglog"))
        afunargs <- list(dep = x$param["dep"])
    if(x$model  %in% c("alog","aneglog"))
        afunargs <- list(dep = x$param["dep"], asy = x$param[c("asy1","asy2")])
    if(x$model  %in% c("bilog","negbilog","ct"))
        afunargs <- list(alpha = x$param["alpha"], beta = x$param["beta"])
    afunargs <- c(list(add = TRUE, lty = lty, model = x$model), afunargs)
    do.call("abvpar", afunargs)
    invisible(x)
}

"mtransform"<- 
function(x, p, inv = FALSE, drp = FALSE)
{
    if(is.list(p)) {
      if(is.null(dim(x)) && length(x) != length(p))
        stop(paste("`p' must have", length(x), "elements"))
      if(!is.null(dim(x)) && ncol(x) != length(p))
        stop(paste("`p' must have", ncol(x), "elements"))
      if(is.null(dim(x))) dim(x) <- c(1, length(p))
      for(i in 1:length(p))
        x[,i] <- Recall(x[,i], p[[i]], inv = inv)
      if(ncol(x) == 1 || (nrow(x) == 1 && drp)) x <- drop(x)
      return(x)
    }
    if(is.null(dim(x))) dim(x) <- c(length(x), 1)
    p <- matrix(t(p), nrow = nrow(x), ncol = 3, byrow = TRUE)
    if(min(p[,2]) <= 0) stop("invalid marginal scale")
    gumind <- (p[,3] == 0)
    nzshapes <- p[!gumind,3]
    if(!inv) {
        x <- (x - p[,1])/p[,2]
        x[gumind, ] <- exp(-x[gumind, ])
        if(any(!gumind))
        x[!gumind, ] <- pmax(1 + nzshapes*x[!gumind, ], 0)^(-1/nzshapes)
    }
    else {
        x[gumind, ] <- p[gumind,1] - p[gumind,2] * log(x[gumind, ])
        x[!gumind, ] <- p[!gumind,1] + p[!gumind,2] *
          (x[!gumind, ]^(-nzshapes) - 1)/nzshapes
    }
    if(ncol(x) == 1 || (nrow(x) == 1 && drp)) x <- drop(x)
    x
}

"bvdepfn" <- 
function(x, a, add, lty, lwd, col, blty, xlab, ylab, xlim, ylim, ...)
{
    if(!add)  { 
        plot(x, a, type="n", xlab = xlab, ylab = ylab,
             xlim = xlim, ylim = ylim, ...) 
        polygon(c(0, 0.5, 1), c(1, 0.5, 1), lty = blty)  
    }
    lines(x, a, lty = lty, lwd = lwd, col = col)
}

na.vals <- function(x) {
    x <- cbind(is.na(x[,1]),2*is.na(x[,2]))
    drop(x %*% c(1,1))
}

"bvstart.vals" <- 
function(x, start, nsloc1, nsloc2, nmdots, param, loc.param1, loc.param2, model, obj = "bvevd", u = NULL)
{
  if(missing(start)) {
    start <- as.list(numeric(length(param)))
    names(start) <- param
    if(obj == "bvevd") {
      st1 <- as.list(fitted(fgev(x[,1], std.err = FALSE, nsloc = nsloc1)))
      st2 <- as.list(fitted(fgev(x[,2], std.err = FALSE, nsloc = nsloc2)))
    }
    if(obj == "bvpot") {
      st1 <- as.list(fitted(fpot(x[,1], u[1], std.err = FALSE)))
      st2 <- as.list(fitted(fpot(x[,2], u[2], std.err = FALSE)))
    }
    start[c(loc.param1, "scale1", "shape1")] <- st1
    tmp2 <- loc.param2
    if("scale2" %in% param) tmp2 <- c(tmp2, "scale2")
    if("shape2" %in% param) tmp2 <- c(tmp2, "shape2")
    tmp <- sub("2", "", tmp2)
    start[tmp2] <- st2[tmp]
    if(model == "log") start[["dep"]] <- 0.75
    if(model == "alog") {
        start[["asy1"]] <- 0.75
        if("asy2" %in% param) start[["asy2"]] <- 0.75 
        start[["dep"]] <- 0.65
    }
    if(model == "hr") start[["dep"]] <- 1
    if(model == "neglog") start[["dep"]] <- 0.6
    if(model == "aneglog") {
        start[["asy1"]] <- 0.75
        if("asy2" %in% param) start[["asy2"]] <- 0.75 
        start[["dep"]] <- 0.8
    }
    if(model == "bilog") start[["alpha"]] <- start[["beta"]] <- 0.75
    if(model == "negbilog") start[["alpha"]] <- start[["beta"]] <- 1/0.6
    if(model == "ct") {
      start[["alpha"]] <- 0.6
      if("beta" %in% param) start[["beta"]] <- 0.6 
    }
    start <- start[!(param %in% nmdots)]
  }
  if(any(!is.na(match(names(start),c("mar1","mar2","asy"))))) {
    if(("mar1" %in% names(start)) && (length(start$mar1) !=
      (2+length(loc.param1)))) stop("mar1 in `start' has incorrect length")
    if(("mar2" %in% names(start)) && (length(start$mar2) !=
      (2+length(loc.param2)))) stop("mar2 in `start' has incorrect length")
    if(("asy" %in% names(start)) && (length(start$asy) != 2))
      stop("asy in `start' should have length two")
    start <- unlist(start)
    names(start)[grep("mar1",names(start))] <- c(loc.param1,"scale1","shape1")
    names(start)[grep("mar2",names(start))] <- c(loc.param2,"scale2","shape2")
    start <- as.list(start)
  }
  if(!is.list(start)) 
    stop("`start' must be a named list")
  start
}

"sep.bvdata" <- 
function(x, obj = "bvevd", u = NULL, censored = TRUE)
{
    if(obj == "bvevd") {
      na <- na.vals(x)
      if(!any(na == 0))
        stop("`x' must have at least one complete observation")
      x.m1 <- as.double(x[na == 2, 1])
      n.m1 <- as.integer(length(x.m1))
      x.m2 <- as.double(x[na == 1, 2])
      n.m2 <- as.integer(length(x.m2))
      x.full <- x[na == 0, , drop = FALSE]
      x1 <- as.double(x.full[,1])
      x2 <- as.double(x.full[,2])
      n <- as.integer(nrow(x.full))
      spx <- list(x.m1 = x.m1, n.m1 = n.m1, x.m2 = x.m2, n.m2 = n.m2,
        x1 = x1, x2 = x2, n = n, na = na)
    }
    if(obj == "bvpot") {
      x1 <- x[,1]
      x2 <- x[,2]
      n <- length(x1)
      r1 <- r2 <- NULL
      iau1 <- (x1 > u[1]) & !is.na(x1)
      iau2 <- (x2 > u[2]) & !is.na(x2)
      nat <- c(sum(iau1), sum(iau2), sum(iau1 & iau2))
      lambda1 <- 1 - sum(!iau1) / n
      lambda2 <- 1 - sum(!iau2) / n
      lambda <- c(lambda1, lambda2)
      if(!censored) {
        x1[is.na(x1)] <- mean(x1[!iau1], na.rm = TRUE)
        x2[is.na(x2)] <- mean(x2[!iau2], na.rm = TRUE)
        r1 <- 1 - rank(x1)/n
        r2 <- 1 - rank(x2)/n
        r1[iau1] <- lambda1
        r2[iau2] <- lambda2
      }
      x1 <- x1 - u[1]
      x2 <- x2 - u[2]
      x1[!iau1] <- 0
      x2[!iau2] <- 0
      i0 <- iau1 | iau2
      x1 <- x1[i0] ; x2 <- x2[i0]
      if(!censored) {
       r1 <- r1[i0] ; r2 <- r2[i0]
      }
      nn <- length(x1)
      thdi <- as.logical(x1) + 2*as.logical(x2)
      spx <- list(x1 = x1, x2 = x2, nn = nn, n = n, thdi = thdi, lambda =
        lambda, r1 = r1, r2 = r2, nat = nat)
    }
    spx
}

"bvpost.optim" <- 
function(x, opt, nm, nsloc1, nsloc2, fixed.param, std.err, dsm, corr, sym, cmar, model)
{
    if (opt$convergence != 0) {
        warning(paste("optimization for", model, "may not have succeeded"))
        if(opt$convergence == 1) opt$convergence <- "iteration limit reached"
    }
    else opt$convergence <- "successful"
    if(std.err) {
        tol <- .Machine$double.eps^0.5
        var.cov <- qr(opt$hessian, tol = tol)
        if(var.cov$rank != ncol(var.cov$qr)) 
            stop(paste("observed information matrix for", model,
                       "is singular; use std.err = FALSE"))
        var.cov <- solve(var.cov, tol = tol)
        std.err <- diag(var.cov)
        if(any(std.err <= 0))
            stop(paste("observed information matrix for", model,
                       "is singular; use std.err = FALSE"))
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
    fixed <- unlist(fixed.param)
    param <- c(opt$par, fixed)
    fixed2 <- NULL
    if(cmar[1]) {
      loc.param1 <- paste("loc1", c("",names(nsloc1)), sep="")
      fixed2 <- c(fixed2, param[loc.param1])
    }
    if(cmar[2]) fixed2 <- c(fixed2, param["scale1"])
    if(cmar[3]) fixed2 <- c(fixed2, param["shape1"])
    if(sym) {
      if(model %in% c("alog","aneglog")) fixed2 <- c(fixed2, param["asy1"])
      if(model == "ct") fixed2 <- c(fixed2, param["alpha"])
    }
    if(!is.null(fixed2)) {
      names(fixed2) <- sub("1", "2", names(fixed2))
      names(fixed2) <- sub("alpha", "beta", names(fixed2))
    }
    param <- c(param, fixed2)
    # Transform to stationarity
    x2 <- x
    if(!is.null(nsloc1)) {
        trend <- param[paste("loc1", names(nsloc1), sep="")]
        trend <- drop(as.matrix(nsloc1) %*% trend)
        x2[,1] <- x[,1] - trend
    }
    if(!is.null(nsloc2)) {
        trend <- param[paste("loc2", names(nsloc2), sep="")]
        trend <- drop(as.matrix(nsloc2) %*% trend)
        x2[,2] <- x[,2] - trend
    }
    # End transform
    # Dependence summary
    if(dsm) {
      dep.sum <- numeric(3)
      if(model %in% c("log", "hr", "neglog")) {
        dep <- param["dep"]
        dep.sum[1] <- 2*(1 - abvpar(dep = dep, model = model))
        dep.sum[2] <- 4 * integrate(function(x) 1-abvpar(x, dep = dep,
          model = model), 0, 1)$value
      }
      if(model %in% c("alog", "aneglog")) {
        dep <- param["dep"]
        asy <- param[c("asy1", "asy2")]
        dep.sum[1] <- 2*(1 - abvpar(dep = dep, asy = asy, model = model))
        dep.sum[2] <- 4 * integrate(function(x) 1-abvpar(x, dep = dep,
          asy = asy, model = model), 0, 1)$value
        dffn <- function(x) abvpar(x, dep = dep, asy = asy, model = model) -
          abvpar(x, dep = dep, asy = rev(asy), model = model)
        dep.sum[3] <- 4*integrate(dffn, 0, .5)$value / (3 - 2*sqrt(2))
      }
      if(model %in% c("bilog", "negbilog", "ct")) {
        alpha <- param["alpha"]
        beta <- param["beta"]
        dep.sum[1] <- 2*(1-abvpar(alpha = alpha, beta = beta, model = model))
        dep.sum[2] <- 4 * integrate(function(x) 1-abvpar(x, alpha = alpha,
          beta = beta, model = model), 0, 1)$value
        dffn <- function(x) abvpar(x, alpha = alpha, beta = beta, model =
          model) - abvpar(x, alpha = beta, beta = alpha, model = model)
        dep.sum[3] <- 4*integrate(dffn, 0, .5)$value / (3 - 2*sqrt(2))
      }
    }
    else dep.sum <- NULL
    # End dependence summary
    list(estimate = opt$par, std.err = std.err, fixed = fixed,
    fixed2 = fixed2, param = param, deviance = 2*opt$value,
    dep.summary = dep.sum, corr = corr, convergence = opt$convergence,
    counts = opt$counts, message = opt$message, data = x, tdata = x2,
    nsloc1 = nsloc1, nsloc2 = nsloc2, n = nrow(x), sym = sym, cmar =
    cmar, model = model)
}






