"rbvlog"<-
# Uses Algorithm 1.1 in Stephenson(2002)
function(n, dep, mar1 = c(0,1,0), mar2 = mar1)
{
    if(length(dep) != 1 || mode(dep) != "numeric" || dep <= 0 ||
        dep > 1) stop("invalid argument for `dep'")
    sim <- .C("rbvlog_shi",
               as.integer(n), as.double(dep), sim = double(2*n))$sim
    bvmtransform(matrix(1/sim, ncol=2, byrow=TRUE), mar1, mar2, inv = TRUE)
}

"rbvalog"<-
# Uses Algorithm 1.2 in Stephenson(2002)
function(n, dep, asy, mar1 = c(0,1,0), mar2 = mar1)
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
              sim = double(2*n))$sim
    bvmtransform(matrix(1/sim, ncol=2, byrow=TRUE), mar1, mar2, inv = TRUE)
}
 
"rbvhr" <-
function(n, dep, mar1 = c(0,1,0), mar2 = mar1)
{
    if(length(dep) != 1 || mode(dep) != "numeric" || dep <= 0)
        stop("invalid argument for `dep'")
    sim <- matrix(runif(2*n),nrow=n,ncol=2)
    condcop <- function(m1,m2,oldm1) {
        tm1 <- -log(m1)
        tm2 <- -log(m2)
        idep <- 1/dep
        v <- tm2 * pnorm(idep + log(tm2/tm1)*dep/2) +
            tm1 * pnorm(idep + log(tm1/tm2)*dep/2)
        exp(-v) * pnorm(idep + log(tm2/tm1)*dep/2) / m2 - oldm1
    }
    for(i in 1:n)
        sim[i,1] <- uniroot(condcop, lower = .Machine$double.eps^0.5,
            upper = 1-.Machine$double.eps^0.5, m2 = sim[i,2],
            oldm1 = sim[i,1])$root
    bvmtransform(-log(sim), mar1, mar2, inv = TRUE)
}

"rbvneglog"<- 
function(n, dep, mar1 = c(0,1,0), mar2 = mar1)
{
    if(length(dep) != 1 || mode(dep) != "numeric" || dep <= 0)
        stop("invalid argument for `dep'")
    sim <- matrix(runif(2*n),nrow=n,ncol=2)
    condcop <- function(m1,m2,oldm1) {
        tm1 <- -log(m1)
        tm2 <- -log(m2)
        idep <- 1/dep
        v <- (tm2^(-dep)+ tm1^(-dep))^(-idep)
        exp(v) * m1 * (1 - (1 + (tm2/tm1)^dep)^(-1-idep)) - oldm1
    }
    for(i in 1:n)
        sim[i,1] <- uniroot(condcop, lower = .Machine$double.eps^0.5,
            upper = 1-.Machine$double.eps^0.5, m2 = sim[i,2],
            oldm1 = sim[i,1])$root
    bvmtransform(-log(sim), mar1, mar2, inv = TRUE)
}

"rbvaneglog"<- 
function(n, dep, asy, mar1 = c(0,1,0), mar2 = mar1)
{
    if(length(dep) != 1 || mode(dep) != "numeric" || dep <= 0)
        stop("invalid argument for `dep'")
    if(length(asy) != 2 || mode(asy) != "numeric" || min(asy) < 0 ||
       max(asy) > 1) stop("invalid argument for `asy'")
    sim <- matrix(runif(2*n),nrow=n,ncol=2)
    condcop <- function(m1,m2,oldm1) {
        tm1 <- -log(m1)
        tm2 <- -log(m2)
        idep <- 1/dep
        v <- ((asy[1]*tm2)^(-dep)+ (asy[2]*tm1)^(-dep))
        exp(v^-idep) * m1 * (1 - asy[1]^(-dep) * tm2^(-dep-1) * v^(-idep-1)) -
        oldm1
    }
    for(i in 1:n)
        sim[i,1] <- uniroot(condcop, lower = .Machine$double.eps^0.5,
            upper = 1-.Machine$double.eps^0.5, m2 = sim[i,2],
            oldm1 = sim[i,1])$root
    bvmtransform(-log(sim), mar1, mar2, inv = TRUE)
}

"pbvlog"<- 
function(q, dep, mar1 = c(0,1,0), mar2 = mar1)
{
    if(length(dep) != 1 || mode(dep) != "numeric" || dep <= 0 ||
        dep > 1) stop("invalid argument for `dep'")
    if(is.null(dim(q))) dim(q) <- c(1,2)
    q <- bvmtransform(q, mar1, mar2)
    v <- apply(q^(1/dep),1,sum)^dep
    exp(-v)
}

"abvlog"<- 
function(x = 0.5, dep, plot = FALSE, border = TRUE, add = FALSE, lty = 1,
         blty = 3, xlim = c(0,1), ylim = c(0.5,1), xlab = "", ylab = "", ...)
{
    if(mode(x) != "numeric" || any(x < 0) || any(x > 1))
        stop("invalid argument for `x'")
    if(length(dep) != 1 || mode(dep) != "numeric" || dep <= 0 ||
        dep > 1) stop("invalid argument for `dep'")
    if(plot || add) x <- seq(0, 1, length = 100)
    idep <- 1/dep
    a <- (x^idep + (1-x)^idep)^dep
    if(plot || add) {
        bvdep(x, a, border, add, lty, blty, xlab, ylab, xlim, ylim, ...)  
        return(invisible(a))
    }
    a
}
  
"pbvalog"<- 
function(q, dep, asy, mar1 = c(0,1,0), mar2 = mar1)
{
    if(length(dep) != 1 || mode(dep) != "numeric" || dep <= 0 ||
        dep > 1) stop("invalid argument for `dep'")
    if(length(asy) != 2 || mode(asy) != "numeric" || min(asy) < 0 ||
       max(asy) > 1) stop("invalid argument for `asy'")
    if(is.null(dim(q))) dim(q) <- c(1,2)
    q <- bvmtransform(q, mar1, mar2)
    asy <- rep(asy,rep(nrow(q),2))
    v <- apply((asy*q)^(1/dep),1,sum)^dep + apply((1-asy)*q,1,sum)
    exp(-v)
}

"abvalog"<- 
function(x = 0.5, dep, asy, plot = FALSE, border = TRUE, add = FALSE, lty = 1,
         blty = 3, xlim = c(0,1), ylim = c(0.5,1), xlab = "", ylab = "", ...)
{
    if(mode(x) != "numeric" || any(x < 0) || any(x > 1))
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
        bvdep(x, a, border, add, lty, blty, xlab, ylab, xlim, ylim, ...) 
        return(invisible(a))
    }
    a
}

"pbvhr" <-
function(q, dep, mar1 = c(0,1,0), mar2 = mar1)
{
    if(length(dep) != 1 || mode(dep) != "numeric" || dep <= 0)
        stop("invalid argument for `dep'")
    if(is.null(dim(q))) dim(q) <- c(1,2)
    q <- bvmtransform(q, mar1, mar2)
    fn <- function(x1,x2) x1*pnorm(1/dep + dep * log(x1/x2) / 2)
    v <- fn(q[,1],q[,2]) + fn(q[,2],q[,1])
    exp(-v)
}

"abvhr" <-
function(x = 0.5, dep, plot = FALSE, border = TRUE, add = FALSE, lty = 1,
         blty = 3, xlim = c(0,1), ylim = c(0.5,1), xlab = "", ylab = "", ...)
{
    if(mode(x) != "numeric" || any(x < 0) || any(x > 1))
        stop("invalid argument for `x'")
    if(length(dep) != 1 || mode(dep) != "numeric" || dep <= 0)
        stop("invalid argument for `dep'")
    if(plot || add) x <- seq(0, 1, length = 100)
    fn <- function(z) z*pnorm(1/dep + dep * log(z/(1-z)) / 2)
    a <- fn(x) + fn(1-x)
    if(plot || add) {
        bvdep(x, a, border, add, lty, blty, xlab, ylab, xlim, ylim, ...)  
        return(invisible(a))
    }
    a
}

"pbvneglog"<- 
function(q, dep, mar1 = c(0,1,0), mar2 = mar1)
{
    if(length(dep) != 1 || mode(dep) != "numeric" || dep <= 0)
        stop("invalid argument for `dep'")
    if(is.null(dim(q))) dim(q) <- c(1,2)
    q <- bvmtransform(q, mar1, mar2)
    v <- apply(q,1,sum) - apply(q^(-dep),1,sum)^(-1/dep)
    exp(-v)
}

"abvneglog"<- 
function(x = 0.5, dep, plot = FALSE, border = TRUE, add = FALSE, lty = 1,
         blty = 3, xlim = c(0,1), ylim = c(0.5,1), xlab = "", ylab = "", ...)
{
    if(mode(x) != "numeric" || any(x < 0) || any(x > 1))
        stop("invalid argument for `x'")  
    if(length(dep) != 1 || mode(dep) != "numeric" || dep <= 0)
        stop("invalid argument for `dep'")
    if(plot || add) x <- seq(0, 1, length = 100)
    a <- 1 - (x^(-dep) + (1-x)^(-dep))^(-1/dep)
    if(plot || add) {
        bvdep(x, a, border, add, lty, blty, xlab, ylab, xlim, ylim, ...) 
        return(invisible(a))
    }
    a
}

"pbvaneglog"<- 
function(q, dep, asy, mar1 = c(0,1,0), mar2 = mar1)
{
    if(length(dep) != 1 || mode(dep) != "numeric" || dep <= 0)
        stop("invalid argument for `dep'")
    if(length(asy) != 2 || mode(asy) != "numeric" || min(asy) < 0 ||
       max(asy) > 1) stop("invalid argument for `asy'")
    if(is.null(dim(q))) dim(q) <- c(1,2)
    q <- bvmtransform(q, mar1, mar2)
    asy <- rep(asy,rep(nrow(q),2))
    v <- apply(q,1,sum) - apply((asy*q)^(-dep),1,sum)^(-1/dep)
    exp(-v)
}

"abvaneglog"<- 
function(x = 0.5, dep, asy, plot = FALSE, border = TRUE, add = FALSE, lty = 1,
         blty = 3, xlim = c(0,1), ylim = c(0.5,1), xlab = "", ylab = "", ...)
{
    if(mode(x) != "numeric" || any(x < 0) || any(x > 1))
        stop("invalid argument for `x'")
    if(length(dep) != 1 || mode(dep) != "numeric" || dep <= 0)
        stop("invalid argument for `dep'")
    if(plot || add) x <- seq(0, 1, length = 100)
    if(length(asy) != 2 || mode(asy) != "numeric" || min(asy) < 0 ||
       max(asy) > 1) stop("invalid argument for `asy'")
    a <- 1 - ((asy[1]*x)^(-dep) + (asy[2]*(1-x))^(-dep))^(-1/dep)
    if(plot || add) {
        bvdep(x, a, border, add, lty, blty, xlab, ylab, xlim, ylim, ...)  
        return(invisible(a))
    }
    a
}

"dbvlog"<- 
function(x, dep, mar1 = c(0,1,0), mar2 = mar1, log = FALSE)
{
    if(length(dep) != 1 || mode(dep) != "numeric" || dep <= 0 ||
        dep > 1) stop("invalid argument for `dep'")
    if(is.null(dim(x))) dim(x) <- c(1,2)
    if(is.null(dim(mar1))) dim(mar1) <- c(1,3)
    if(is.null(dim(mar2))) dim(mar2) <- c(1,3)
    d <- numeric(nrow(x))
    x <- bvmtransform(x, mar1, mar2)
    ext <- apply(x,1,function(z) any(z %in% c(0,Inf)))
    d[ext] <- -Inf
    if(sum(!ext)) {
        x <- x[!ext, ,drop=FALSE]
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
function(x, dep, asy, mar1 = c(0,1,0), mar2 = mar1, log = FALSE)
{
    if(length(dep) != 1 || mode(dep) != "numeric" || dep <= 0 ||
        dep > 1) stop("invalid argument for `dep'")
    if(length(asy) != 2 || mode(asy) != "numeric" || min(asy) < 0 ||
       max(asy) > 1) stop("invalid argument for `asy'")
    if(is.null(dim(x))) dim(x) <- c(1,2)
    if(is.null(dim(mar1))) dim(mar1) <- c(1,3)
    if(is.null(dim(mar2))) dim(mar2) <- c(1,3)
    d <- numeric(nrow(x))
    x <- bvmtransform(x, mar1, mar2)
    ext <- apply(x,1,function(z) any(z %in% c(0,Inf)))
    d[ext] <- -Inf
    if(sum(!ext)) {
        x <- x[!ext, ,drop=FALSE]
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
    if(is.null(dim(mar1))) dim(mar1) <- c(1,3)
    if(is.null(dim(mar2))) dim(mar2) <- c(1,3)
    d <- numeric(nrow(x))
    x <- bvmtransform(x, mar1, mar2)
    ext <- apply(x,1,function(z) any(z %in% c(0,Inf)))
    d[ext] <- -Inf
    if(sum(!ext)) {
        x <- x[!ext, ,drop=FALSE]
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
    if(is.null(dim(mar1))) dim(mar1) <- c(1,3)
    if(is.null(dim(mar2))) dim(mar2) <- c(1,3)
    d <- numeric(nrow(x))
    x <- bvmtransform(x, mar1, mar2)
    ext <- apply(x,1,function(z) any(z %in% c(0,Inf)))
    d[ext] <- -Inf
    if(sum(!ext)) {
        x <- x[!ext, ,drop=FALSE]
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
function(x, dep, asy, mar1 = c(0,1,0), mar2 = mar1, log = FALSE)
{
    if(length(dep) != 1 || mode(dep) != "numeric" || dep <= 0)
        stop("invalid argument for `dep'")
    if(length(asy) != 2 || mode(asy) != "numeric" || min(asy) < 0 ||
       max(asy) > 1) stop("invalid argument for `asy'")
    if(is.null(dim(x))) dim(x) <- c(1,2)
    if(is.null(dim(mar1))) dim(mar1) <- c(1,3)
    if(is.null(dim(mar2))) dim(mar2) <- c(1,3)
    d <- numeric(nrow(x))
    x <- bvmtransform(x, mar1, mar2)
    ext <- apply(x,1,function(z) any(z %in% c(0,Inf)))
    d[ext] <- -Inf
    if(sum(!ext)) {
        x <- x[!ext, ,drop=FALSE]
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

"fbvlog"<- 
function(x, start, ...)
{
    if(any(!is.na(match(names(start),c("mar1","mar2"))))) {
        start <- unlist(start)
        names(start)[grep("mar1",names(start))] <- c("loc1","scale1","shape1")
        names(start)[grep("mar2",names(start))] <- c("loc2","scale2","shape2")
        start <- as.list(start)
    }
    ddbvlog <- function(x, dep, loc1 = 0, scale1 = 1, shape1 = 0,
                    loc2 = loc1, scale2 = scale1, shape2 = shape1, log) {
                if(any(c(dep,scale1,scale2) <= 0) || dep > 1)
                    return(-Inf)
                dbvlog(x = x, dep = dep,
                mar1 = c(loc1,scale1,shape1), mar2 = c(loc2,scale2,shape2),
                log = log)
                }
    fdensfun(x = x, densfun = ddbvlog, start = start, ...)
}

"fbvalog"<- 
function(x, start, ...)
{
    if(any(!is.na(match(names(start),c("asy","mar1","mar2"))))) {
        start <- unlist(start)
        names(start)[grep("mar1",names(start))] <- c("loc1","scale1","shape1")
        names(start)[grep("mar2",names(start))] <- c("loc2","scale2","shape2")
        start <- as.list(start)
    }
    ddbvalog <- function(x, dep, asy1, asy2, loc1 = 0, scale1 = 1, shape1 = 0,
                    loc2 = loc1, scale2 = scale1, shape2 = shape1, log) {
                if(any(c(dep,scale1,scale2) <= 0) ||
                   any(c(dep,asy1,asy2) > 1) || any(c(asy1,asy2) < 0))
                    return(-Inf)
                dbvalog(x = x, dep = dep, asy = c(asy1,asy2),
                mar1 = c(loc1,scale1,shape1), mar2 = c(loc2,scale2,shape2),
                log = log)
                }
    fdensfun(x = x, densfun = ddbvalog, start = start, ...)
}

"fbvhr"<- 
function(x, start, ...)
{
    if(any(!is.na(match(names(start),c("mar1","mar2"))))) {
        start <- unlist(start)
        names(start)[grep("mar1",names(start))] <- c("loc1","scale1","shape1")
        names(start)[grep("mar2",names(start))] <- c("loc2","scale2","shape2")
        start <- as.list(start)
    }
    ddbvhr <- function(x, dep, loc1 = 0, scale1 = 1, shape1 = 0,
                  loc2 = loc1, scale2 = scale1, shape2 = shape1, log) {
              if(any(c(dep,scale1,scale2) <= 0))
                  return(-Inf)
              dbvhr(x = x, dep = dep,
              mar1 = c(loc1,scale1,shape1), mar2 = c(loc2,scale2,shape2),
              log = log)
              }
    fdensfun(x = x, densfun = ddbvhr, start = start, ...)
}

"fbvneglog"<- 
function(x, start, ...)
{
    if(any(!is.na(match(names(start),c("mar1","mar2"))))) {
        start <- unlist(start)
        names(start)[grep("mar1",names(start))] <- c("loc1","scale1","shape1")
        names(start)[grep("mar2",names(start))] <- c("loc2","scale2","shape2")
        start <- as.list(start)
    }
    ddbvneglog <- function(x, dep, loc1 = 0, scale1 = 1, shape1 = 0,
                    loc2 = loc1, scale2 = scale1, shape2 = shape1, log) {
                if(any(c(dep,scale1,scale2) <= 0))
                    return(-Inf)
                dbvneglog(x = x, dep = dep,
                mar1 = c(loc1,scale1,shape1), mar2 = c(loc2,scale2,shape2),
                log = log)
                }
    fdensfun(x = x, densfun = ddbvneglog, start = start, ...)
}

"fbvaneglog"<- 
function(x, start, ...)
{
    if(any(!is.na(match(names(start),c("asy","mar1","mar2"))))) {
        start <- unlist(start)
        names(start)[grep("mar1",names(start))] <- c("loc1","scale1","shape1")
        names(start)[grep("mar2",names(start))] <- c("loc2","scale2","shape2")
        start <- as.list(start)
    }
    ddbvaneglog <- function(x, dep, asy1, asy2, loc1 = 0, scale1 = 1,
                       shape1 = 0, loc2 = loc1, scale2 = scale1,
                       shape2 = shape1, log) {
        if(any(c(dep,scale1,scale2) <= 0) || any(c(asy1,asy2) > 1) ||
           any(c(asy1,asy2) < 0))
            return(-Inf)
        dbvaneglog(x = x, dep = dep, asy = c(asy1,asy2),
            mar1 = c(loc1,scale1,shape1), mar2 = c(loc2,scale2,shape2),
            log = log)
                }
    fdensfun(x = x, densfun = ddbvaneglog, start = start, ...)
}

"mtransform"<- 
function(x, p, inv = FALSE)
{
    if(is.null(dim(x))) dim(x) <- c(length(x),1)
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
    if(ncol(x) == 1) dim(x) <- NULL
    x
}

"bvmtransform"<- 
function(x, p1, p2, inv = FALSE)
{
    if(is.null(dim(x))) dim(x) <- c(1,2)
    x[,1] <- mtransform(x[,1], p1, inv = inv)
    x[,2] <- mtransform(x[,2], p2, inv = inv)
    x
}

"bvdep" <- 
function(x, a, border, add, lty, blty, xlab, ylab, xlim, ylim, ...)
{
    if(!add)  { 
        oldpty <- par(pty="s") 
        on.exit(par(oldpty)) 
        plot(x, a, type="n", xlab = xlab, ylab = ylab,
             xlim = xlim, ylim = ylim, ...) 
        if(border)  polygon(c(0, 0.5, 1), c(1, 0.5, 1), lty = blty)  
    }
    lines(x, a, lty = lty)
}


