"rbvlog"<-
# Uses Algorithm 1.1 in Stephenson(2002)
function(n, dep, mar1 = c(0,1,0), mar2 = mar1)
{
    if(length(dep) != 1 || mode(dep) != "numeric" || dep <= 0 ||
        dep > 1) stop("invalid argument for `dep'")
    sim <- .C("rbvlog_shi",
               as.integer(n), as.double(dep), sim = double(2*n),
               PACKAGE = "evd")$sim
    bvmtransform(matrix(1/sim, ncol=2, byrow=TRUE), mar1, mar2, inv = TRUE)
}

"rbvalog"<-
# Uses Algorithm 1.2 in Stephenson(2002)
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
    bvmtransform(matrix(1/sim, ncol=2, byrow=TRUE), mar1, mar2, inv = TRUE)
}

"rbvhr" <-
function(n, dep, mar1 = c(0,1,0), mar2 = mar1)
{
    if(length(dep) != 1 || mode(dep) != "numeric" || dep <= 0)
        stop("invalid argument for `dep'")
    sim <- .C("rbvhr",
        as.double(n), as.double(dep), sim = runif(2*n),
        PACKAGE = "evd")$sim
    sim <- matrix(sim, nrow = n, ncol = 2, byrow = TRUE)
    bvmtransform(-log(sim), mar1, mar2, inv = TRUE)
}

"rbvneglog"<- 
function(n, dep, mar1 = c(0,1,0), mar2 = mar1)
{
    if(length(dep) != 1 || mode(dep) != "numeric" || dep <= 0)
        stop("invalid argument for `dep'")
    sim <- .C("rbvneglog",
        as.double(n), as.double(dep), sim = runif(2*n),
        PACKAGE = "evd")$sim
    sim <- matrix(sim, nrow = n, ncol = 2, byrow = TRUE)
    bvmtransform(-log(sim), mar1, mar2, inv = TRUE)
}

"rbvaneglog"<- 
function(n, dep, asy = c(1,1), mar1 = c(0,1,0), mar2 = mar1)
{
    if(length(dep) != 1 || mode(dep) != "numeric" || dep <= 0)
        stop("invalid argument for `dep'")
    if(length(asy) != 2 || mode(asy) != "numeric" || min(asy) < 0 ||
       max(asy) > 1) stop("invalid argument for `asy'")
    sim <- .C("rbvaneglog",
        as.double(n), as.double(dep), as.double(asy), sim = runif(2*n),
        PACKAGE = "evd")$sim
    sim <- matrix(sim, nrow = n, ncol = 2, byrow = TRUE)
    bvmtransform(-log(sim), mar1, mar2, inv = TRUE)
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
        as.double(n), as.double(alpha), as.double(beta), sim = runif(2*n),
        PACKAGE = "evd")$sim
    sim <- matrix(sim, nrow = n, ncol = 2, byrow = TRUE)
    bvmtransform(-log(sim), mar1, mar2, inv = TRUE)
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
        as.double(n), as.double(alpha), as.double(beta), sim = runif(2*n),
        PACKAGE = "evd")$sim
    sim <- matrix(sim, nrow = n, ncol = 2, byrow = TRUE)
    bvmtransform(-log(sim), mar1, mar2, inv = TRUE)
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
        as.double(n), as.double(alpha), as.double(beta), sim = runif(2*n),
        PACKAGE = "evd")$sim
    sim <- matrix(sim, nrow = n, ncol = 2, byrow = TRUE)
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

"pbvalog"<- 
function(q, dep, asy = c(1,1), mar1 = c(0,1,0), mar2 = mar1)
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

"pbvaneglog"<- 
function(q, dep, asy = c(1,1), mar1 = c(0,1,0), mar2 = mar1)
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

"pbvbilog"<- 
function(q, alpha, beta, mar1 = c(0,1,0), mar2 = mar1)
{
    if(length(alpha) != 1 || mode(alpha) != "numeric")
        stop("invalid argument for `alpha'")
    if(length(beta) != 1 || mode(beta) != "numeric")
        stop("invalid argument for `beta'")
    if(any(c(alpha,beta) <= 0) || any(c(alpha,beta) >= 1))
        stop("`alpha' and `beta' must be in the open interval (0,1)")
    if(is.null(dim(q))) dim(q) <- c(1,2)
    q <- bvmtransform(q, mar1, mar2)
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
    exp(-v)
}

"pbvnegbilog"<- 
function(q, alpha, beta, mar1 = c(0,1,0), mar2 = mar1)
{
    if(length(alpha) != 1 || mode(alpha) != "numeric")
        stop("invalid argument for `alpha'")
    if(length(beta) != 1 || mode(beta) != "numeric")
        stop("invalid argument for `beta'")
    if(any(c(alpha,beta) <= 0))
        stop("`alpha' and `beta' must be non-negative")
    if(is.null(dim(q))) dim(q) <- c(1,2)
    q <- bvmtransform(q, mar1, mar2)
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
    exp(-v)
}

"pbvct" <-
function(q, alpha, beta, mar1 = c(0,1,0), mar2 = mar1)
{
    if(length(alpha) != 1 || mode(alpha) != "numeric")
        stop("invalid argument for `alpha'")
    if(length(beta) != 1 || mode(beta) != "numeric")
        stop("invalid argument for `beta'")
    if(any(c(alpha,beta) <= 0))
        stop("`alpha' and `beta' must be non-negative")
    if(is.null(dim(q))) dim(q) <- c(1,2)
    q <- bvmtransform(q, mar1, mar2)
    print(q)
    u <- (alpha * q[,2]) / (alpha * q[,2] + beta * q[,1])  
    v <- q[,2] * pbeta(u, shape1 = alpha, shape2 = beta + 1) +
      q[,1] * pbeta(u, shape1 = alpha + 1, shape2 = beta, lower.tail = FALSE)
    v[is.infinite(q[,1]) || is.infinite(q[,2])] <- Inf
    v[(q[,1] + q[,2]) == 0] <- 0
    exp(-v)
}

"abvlog"<- 
function(x = 0.5, dep, plot = FALSE, border = TRUE, add = FALSE, lty = 1,
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
        bvdep(x, a, border, add, lty, blty, xlab, ylab, xlim, ylim, ...)  
        return(invisible(a))
    }
    a
}

"abvalog"<- 
function(x = 0.5, dep, asy = c(1,1), plot = FALSE, border = TRUE, add = FALSE,
         lty = 1, blty = 3, xlim = c(0,1), ylim = c(0.5,1), xlab = "",
         ylab = "", ...)
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
        bvdep(x, a, border, add, lty, blty, xlab, ylab, xlim, ylim, ...) 
        return(invisible(a))
    }
    a
}

"abvhr" <-
function(x = 0.5, dep, plot = FALSE, border = TRUE, add = FALSE, lty = 1,
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
        bvdep(x, a, border, add, lty, blty, xlab, ylab, xlim, ylim, ...)  
        return(invisible(a))
    }
    a
}

"abvneglog"<- 
function(x = 0.5, dep, plot = FALSE, border = TRUE, add = FALSE, lty = 1,
         blty = 3, xlim = c(0,1), ylim = c(0.5,1), xlab = "", ylab = "", ...)
{
    if(mode(x) != "numeric" || any(x < 0,na.rm=TRUE) || any(x > 1,na.rm=TRUE))
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

"abvaneglog"<- 
function(x = 0.5, dep, asy = c(1,1), plot = FALSE, border = TRUE, add = FALSE,
         lty = 1, blty = 3, xlim = c(0,1), ylim = c(0.5,1), xlab = "",
         ylab = "", ...)
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
        bvdep(x, a, border, add, lty, blty, xlab, ylab, xlim, ylim, ...)  
        return(invisible(a))
    }
    a
}

"abvbilog"<- 
function(x = 0.5, alpha, beta, plot = FALSE, border = TRUE, add = FALSE,
         lty = 1, blty = 3, xlim = c(0,1), ylim = c(0.5,1), xlab = "",
         ylab = "", ...)
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
        bvdep(x, a, border, add, lty, blty, xlab, ylab, xlim, ylim, ...)  
        return(invisible(a))
    }
    a
}

"abvnegbilog"<- 
function(x = 0.5, alpha, beta, plot = FALSE, border = TRUE, add = FALSE,
         lty = 1, blty = 3, xlim = c(0,1), ylim = c(0.5,1), xlab = "",
         ylab = "", ...)
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
        bvdep(x, a, border, add, lty, blty, xlab, ylab, xlim, ylim, ...)  
        return(invisible(a))
    }
    a
}

"abvct" <-
function(x = 0.5, alpha, beta, plot = FALSE, border = TRUE, add = FALSE,
         lty = 1, blty = 3, xlim = c(0,1), ylim = c(0.5,1), xlab = "",
         ylab = "", ...)
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
        bvdep(x, a, border, add, lty, blty, xlab, ylab, xlim, ylim, ...)  
        return(invisible(a))
    }
    a
}

"abvnonpar"<- 
function(x = 0.5, data, nsloc1 = NULL, nsloc2 = NULL,
         method = c("cfg","deheuvels","pickands"), modify = 0,
         wf = function(t) t, plot = FALSE, border = TRUE, add = FALSE,
         lty = 1, blty = 3, xlim = c(0,1), ylim = c(0.5,1),
         xlab = "", ylab = "", ...)
{
    if(mode(x) != "numeric" || any(x < 0,na.rm=TRUE) || any(x > 1,na.rm=TRUE))
        stop("invalid argument for `x'")
    if(!is.null(nsloc1)) {
        if(is.vector(nsloc1))
            nsloc1 <- data.frame(trend = nsloc1)
        if(!is.data.frame(nsloc1))
            stop("`nsloc1' must be a vector or data frame")
        if(nrow(nsloc1) != nrow(data))
            stop("`nsloc1' and `data' are not compatible")
        nslocmat1 <- cbind(1,as.matrix(nsloc1))
    }
    if(!is.null(nsloc2)) {
        if(is.vector(nsloc2))
            nsloc2 <- data.frame(trend = nsloc2)
        if(!is.data.frame(nsloc2))
            stop("`nsloc2' must be a vector or data frame")
        if(nrow(nsloc2) != nrow(data))
            stop("`nsloc2' and `data' are not compatible")
        nslocmat2 <- cbind(1,as.matrix(nsloc2))
    }
    mle.m1 <- fgev(data[,1], nsloc = nsloc1, std.err = FALSE)$estimate
    loc.mle.m1 <- mle.m1[grep("^loc", names(mle.m1))]
    if(is.null(nsloc1))
        loc.mle.m1 <- rep(loc.mle.m1, nrow(data))
    else loc.mle.m1 <- nslocmat1 %*% loc.mle.m1
    mle.m1 <- cbind(loc.mle.m1, mle.m1["scale"], mle.m1["shape"])
    mle.m2 <- fgev(data[,2], nsloc = nsloc2, std.err = FALSE)$estimate
    loc.mle.m2 <- mle.m2[grep("^loc", names(mle.m2))]
    if(is.null(nsloc2))
        loc.mle.m2 <- rep(loc.mle.m2, nrow(data))
    else loc.mle.m2 <- nslocmat2 %*% loc.mle.m2
    mle.m2 <- cbind(loc.mle.m2, mle.m2["scale"], mle.m2["shape"])
    data <- bvmtransform(data, mle.m1, mle.m2)
    data <- na.omit(data)
    
    if(plot || add) x <- seq(0, 1, length = 100)
    method <- match.arg(method)
    nn <- nrow(data)
    
    if(method == "pickands") {
        if(modify != 2) {
            a <- numeric(length(x))
            for(i in 1:length(x))
                a[i] <- sum(pmin(data[,1]/x[i], data[,2]/(1-x[i])))
            a <- nn / a
        }
        if(modify == 1) a <- pmin(1, pmax(a, x, 1-x))
        if(modify == 2) {
            x2 <- seq(0, 1, length = 250)
            a <- numeric(250)
            for(i in 1:250)
                a[i] <- sum(pmin(data[,1]/x2[i], data[,2]/(1-x2[i])))
            a <- nn / a
            a <- pmin(1, pmax(a, x2, 1-x2))
            inch <- chull(x2, a)
            a <- a[inch] ; x2 <- x2[inch]
            a <- approx(x2, a, xout = x, method="linear")$y
        }
    } 
    if(method == "deheuvels") {
        if(modify != 2) {
            a <- numeric(length(x))
            for(i in 1:length(x))
                a[i] <- sum(pmin(data[,1]/x[i], data[,2]/(1-x[i])))
            a <- nn / (a - x * sum(data[,1]) - (1-x) * sum(data[,2]) + nn)
        }
        if(modify == 1) a <- pmin(1, pmax(a, x, 1-x))
        if(modify == 2) {
            x2 <- seq(0, 1, length = 250)
            a <- numeric(250)
            for(i in 1:250)
                a[i] <- sum(pmin(data[,1]/x2[i], data[,2]/(1-x2[i]))) 
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
        if(modify != 2) {
            step.pos <- as.numeric(cut(1-x, breaks = c(-0.1, zvec, 1))) - 1
            step.pos.r <- step.pos / nn
            a <- (1-x)^step.pos.r * x^(1 - step.pos.r) *
                qvec[nn]^(1 - wf(x)) / ifelse(step.pos, qvec[step.pos], 1)
        }
        if(modify == 1) a <- pmin(1, pmax(a, x, 1-x))
        if(modify == 2) {
            x2 <- seq(0, 1, length = 250)
            step.pos <- as.numeric(cut(1-x2, breaks = c(-0.1, zvec, 1))) - 1
            step.pos.r <- step.pos / nn
            a <- (1-x2)^step.pos.r * x2^(1 - step.pos.r) *
                qvec[nn]^(1 - wf(x2)) / ifelse(step.pos, qvec[step.pos], 1)
            a <- pmin(1, pmax(a, x2, 1-x2))
            inch <- chull(x2, a)
            a <- a[inch] ; x2 <- x2[inch]
            a <- approx(x2, a, xout = x, method="linear")$y
        }
    }
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
    mar1 <- matrix(t(mar1), nrow = nrow(x), ncol = 3, byrow = TRUE)
    mar2 <- matrix(t(mar2), nrow = nrow(x), ncol = 3, byrow = TRUE)
    d <- numeric(nrow(x))
    x <- bvmtransform(x, mar1, mar2)
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
    x <- bvmtransform(x, mar1, mar2)
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
    x <- bvmtransform(x, mar1, mar2)
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
    x <- bvmtransform(x, mar1, mar2)
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
    x <- bvmtransform(x, mar1, mar2)
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
    x <- bvmtransform(x, mar1, mar2)
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
    x <- bvmtransform(x, mar1, mar2)
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
    x <- bvmtransform(x, mar1, mar2)
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
    
"fbvlog"<- 
function(x, start, ..., nsloc1 = NULL, nsloc2 = NULL, std.err = TRUE,
         method = "BFGS")
{
    nlbvlog <- function(loc1 = 0, scale1 = 1, shape1 = 0, loc2 = loc1,
                        scale2 = scale1, shape2 = shape1, dep)
    {
        if(any(c(scale1,scale2) < 0.01) || dep < 0.05 || dep > 1)
            return(1e6)
        if(!is.null(nsloc1)) {
            ns <- numeric(length(loc.param1))
            for(i in 1:length(ns))
                ns[i] <- get(loc.param1[i])
            loc1 <- drop(nslocmat1 %*% ns)
        }
        else loc1 <- rep(loc1, length.out = length(naind))
        if(!is.null(nsloc2)) {
            ns <- numeric(length(loc.param2))
            for(i in 1:length(ns))
                ns[i] <- get(loc.param2[i])
            loc2 <- drop(nslocmat2 %*% ns)
        }
        else loc2 <- rep(loc2, length.out = length(naind))
        if(any(naind == 2))
            m1l <- .C("nlgev",
                x.m1, n.m1,
                as.double(loc1[naind == 2]), as.double(scale1),
                as.double(shape1), dns = double(1), PACKAGE = "evd")$dns
        else m1l <- 0
        if(is.nan(m1l))
            stop("NaN returned in likelihood evaluation\n\tPlease contact me at a.stephenson@lancaster.ac.uk")
        if(any(naind == 1))
            m2l <- .C("nlgev",
                x.m2, n.m2,
                as.double(loc2[naind == 1]), as.double(scale2),
                as.double(shape2), dns = double(1), PACKAGE = "evd")$dns
        else m2l <- 0
        if(is.nan(m2l))
            stop("NaN returned in likelihood evaluation\n\tPlease contact me at a.stephenson@lancaster.ac.uk")
        bvl <- .C("nlbvlog",
                x1, x2, n,
                as.double(dep), as.double(loc1[naind == 0]), as.double(scale1),
                as.double(shape1), as.double(loc2[naind == 0]),
                as.double(scale2), as.double(shape2), dns = double(1),
                PACKAGE = "evd")$dns
        if(is.nan(bvl))
            stop("NaN returned in likelihood evaluation\n\tPlease contact me at a.stephenson@lancaster.ac.uk")
        if(any(c(m1l,m2l,bvl) == 1e6)) return(1e6)
        else return(m1l + m2l + bvl)
    }
    naind <- na.vals(x)
    if(!any(naind == 0))
        stop("`x' must have at least one complete observation")
    if(!is.null(nsloc1)) {
        if(is.vector(nsloc1))
            nsloc1 <- data.frame(trend = nsloc1)
        if(!is.data.frame(nsloc1))
            stop("`nsloc1' must be a vector or data frame")
        if(nrow(nsloc1) != nrow(x))
            stop("`nsloc1' and `x' are not compatible")
        nslocmat1 <- cbind(1,as.matrix(nsloc1))
    }
    if(!is.null(nsloc2)) {
        if(is.vector(nsloc2))
            nsloc2 <- data.frame(trend = nsloc2)
        if(!is.data.frame(nsloc2))
            stop("`nsloc2' must be a vector or data frame")
        if(nrow(nsloc2) != nrow(x))
            stop("`nsloc2' and `x' are not compatible")
        nslocmat2 <- cbind(1,as.matrix(nsloc2))
    }
    loc.param1 <- paste("loc1", c("",names(nsloc1)), sep="")
    loc.param2 <- paste("loc2", c("",names(nsloc2)), sep="")
    param <- c(loc.param1, "scale1", "shape1", loc.param2, "scale2",
               "shape2", "dep")
    if(missing(start)) {
        start <- as.list(numeric(length(param)))
        names(start) <- param
        st1 <- as.list(fgev(x[,1], std.err = FALSE, nsloc = nsloc1)$estimate)
        st2 <- as.list(fgev(x[,2], std.err = FALSE, nsloc = nsloc2)$estimate)
        start[c(loc.param1, "scale1", "shape1")] <- st1
        start[c(loc.param2, "scale2", "shape2")] <- st2
        start[["dep"]] <- 0.75
        start <- start[!(param %in% names(list(...)))]
    }
    if(any(!is.na(match(names(start),c("mar1","mar2"))))) {
        start <- unlist(start)
        names(start)[grep("mar1",names(start))] <- c("loc1","scale1","shape1")
        names(start)[grep("mar2",names(start))] <- c("loc2","scale2","shape2")
        start <- as.list(start)
    }
    if(!is.list(start)) 
        stop("`start' must be a named list")
    x.m1 <- as.double(x[naind == 2, 1])
    n.m1 <- as.integer(length(x.m1))
    x.m2 <- as.double(x[naind == 1, 2])
    n.m2 <- as.integer(length(x.m2))
    x <- x[naind == 0, , drop = FALSE]
    x1 <- as.double(x[,1])
    x2 <- as.double(x[,2])
    n <- as.integer(nrow(x))
    nm <- names(start)
    l <- length(nm)
    f <- c(as.list(numeric(length(loc.param1))), formals(nlbvlog)[2:3],
           as.list(numeric(length(loc.param2))), formals(nlbvlog)[5:7])
    names(f) <- param
    m <- match(nm, param)
    if(any(is.na(m))) 
        stop("`start' specifies unknown arguments")    
    formals(nlbvlog) <- c(f[m], f[-m])
    nllh <- function(p, ...) nlbvlog(p, ...)
    if(l > 1)
        body(nllh) <- parse(text = paste("nlbvlog(", paste("p[",1:l,
            "]", collapse = ", "), ", ...)"))
    fixed.param <- list(...)[names(list(...)) %in% param]
    start.arg <- c(list(p = unlist(start)), fixed.param)
    if(do.call("nllh", start.arg) == 1e6)
        warning("negative log-likelihood is infinite at starting values")
    
    opt <- optim(start, nllh, hessian = TRUE, ..., method = method)
    if (opt$convergence > 0) 
        warning("optimization may not have succeeded")
    if(std.err) {
        tol <- .Machine$double.eps^0.5
        std.err <- qr(opt$hessian, tol = tol)
        if (std.err$rank != ncol(std.err$qr)) 
            stop("observed information matrix is singular; use std.err = FALSE")
        std.err <- diag(solve(std.err, tol = tol))
        if(any(std.err <= 0))
            stop("observed information matrix is singular; use std.err = FALSE")
        std.err <- sqrt(std.err)
    }
    else std.err <- rep(NA,l)
    names(std.err) <- nm
    list(estimate = opt$par, std.err = std.err, deviance = 2*opt$value,
         counts = opt$counts)
}

"fbvalog"<- 
function(x, start, ..., nsloc1 = NULL, nsloc2 = NULL, std.err = TRUE,
         method = "BFGS")
{
    nlbvalog <- function(loc1 = 0, scale1 = 1, shape1 = 0, loc2 = loc1,
                scale2 = scale1, shape2 = shape1, asy1 = 1, asy2 = 1, dep)
    {
        if(any(c(scale1,scale2) < 0.01) || any(c(dep,asy1,asy2) > 1) ||
           any(c(asy1,asy2) < 0) || dep < 0.05)
            return(1e6)
        if(!is.null(nsloc1)) {
            ns <- numeric(length(loc.param1))
            for(i in 1:length(ns))
                ns[i] <- get(loc.param1[i])
            loc1 <- drop(nslocmat1 %*% ns)
        }
        else loc1 <- rep(loc1, length.out = length(naind))
        if(!is.null(nsloc2)) {
            ns <- numeric(length(loc.param2))
            for(i in 1:length(ns))
                ns[i] <- get(loc.param2[i])
            loc2 <- drop(nslocmat2 %*% ns)
        }
        else loc2 <- rep(loc2, length.out = length(naind))
        if(any(naind == 2))
            m1l <- .C("nlgev",
                x.m1, n.m1,
                as.double(loc1[naind == 2]), as.double(scale1),
                as.double(shape1), dns = double(1), PACKAGE = "evd")$dns
        else m1l <- 0
        if(is.nan(m1l))
            stop("NaN returned in likelihood evaluation\n\tPlease contact me at a.stephenson@lancaster.ac.uk")
        if(any(naind == 1))
            m2l <- .C("nlgev",
                x.m2, n.m2,
                as.double(loc2[naind == 1]), as.double(scale2),
                as.double(shape2), dns = double(1), PACKAGE = "evd")$dns
        else m2l <- 0
        if(is.nan(m2l))
            stop("NaN returned in likelihood evaluation\n\tPlease contact me at a.stephenson@lancaster.ac.uk")
        bvl <- .C("nlbvalog",
                x1, x2, n,
                as.double(dep), as.double(asy1), as.double(asy2),
                as.double(loc1[naind == 0]), as.double(scale1),
                as.double(shape1), as.double(loc2[naind == 0]),
                as.double(scale2), as.double(shape2), dns = double(1),
                PACKAGE = "evd")$dns
        if(is.nan(bvl))
            stop("NaN returned in likelihood evaluation\n\tPlease contact me at a.stephenson@lancaster.ac.uk")
        if(any(c(m1l,m2l,bvl) == 1e6)) return(1e6)
        else return(m1l + m2l + bvl)
    }
    naind <- na.vals(x)
    if(!any(naind == 0))
        stop("`x' must have at least one complete observation")
    if(!is.null(nsloc1)) {
        if(is.vector(nsloc1))
            nsloc1 <- data.frame(trend = nsloc1)
        if(!is.data.frame(nsloc1))
            stop("`nsloc1' must be a vector or data frame")
        if(nrow(nsloc1) != nrow(x))
            stop("`nsloc1' and `x' are not compatible")
        nslocmat1 <- cbind(1,as.matrix(nsloc1))
    }
    if(!is.null(nsloc2)) {
        if(is.vector(nsloc2))
            nsloc2 <- data.frame(trend = nsloc2)
        if(!is.data.frame(nsloc2))
            stop("`nsloc2' must be a vector or data frame")
        if(nrow(nsloc2) != nrow(x))
            stop("`nsloc2' and `x' are not compatible")
        nslocmat2 <- cbind(1,as.matrix(nsloc2))
    }
    loc.param1 <- paste("loc1", c("",names(nsloc1)), sep="")
    loc.param2 <- paste("loc2", c("",names(nsloc2)), sep="")
    param <- c(loc.param1, "scale1", "shape1", loc.param2, "scale2",
               "shape2", "asy1", "asy2", "dep")
    if(missing(start)) {
        start <- as.list(numeric(length(param)))
        names(start) <- param
        st1 <- as.list(fgev(x[,1], std.err = FALSE, nsloc = nsloc1)$estimate)
        st2 <- as.list(fgev(x[,2], std.err = FALSE, nsloc = nsloc2)$estimate)
        start[c(loc.param1, "scale1", "shape1")] <- st1
        start[c(loc.param2, "scale2", "shape2")] <- st2
        start[["asy1"]] <- 0.75
        start[["asy2"]] <- 0.75
        start[["dep"]] <- 0.65
        start <- start[!(param %in% names(list(...)))]
    }
    if(any(!is.na(match(names(start),c("asy","mar1","mar2"))))) {
        start <- unlist(start)
        names(start)[grep("mar1",names(start))] <- c("loc1","scale1","shape1")
        names(start)[grep("mar2",names(start))] <- c("loc2","scale2","shape2")
        start <- as.list(start)
    }
    if(!is.list(start)) 
        stop("`start' must be a named list")
    x.m1 <- as.double(x[naind == 2, 1])
    n.m1 <- as.integer(length(x.m1))
    x.m2 <- as.double(x[naind == 1, 2])
    n.m2 <- as.integer(length(x.m2))
    x <- x[naind == 0, , drop = FALSE]
    x1 <- as.double(x[,1])
    x2 <- as.double(x[,2])
    n <- as.integer(nrow(x))
    nm <- names(start)
    l <- length(nm)
    f <- c(as.list(numeric(length(loc.param1))), formals(nlbvalog)[2:3],
           as.list(numeric(length(loc.param2))), formals(nlbvalog)[5:9])
    names(f) <- param
    
    m <- match(nm, param)
    if(any(is.na(m))) 
        stop("`start' specifies unknown arguments")    
    formals(nlbvalog) <- c(f[m], f[-m])
    nllh <- function(p, ...) nlbvalog(p, ...)
    if(l > 1)
        body(nllh) <- parse(text = paste("nlbvalog(", paste("p[",1:l,
            "]", collapse = ", "), ", ...)"))
    fixed.param <- list(...)[names(list(...)) %in% param]
    start.arg <- c(list(p = unlist(start)), fixed.param)
    if(do.call("nllh", start.arg) == 1e6)
        warning("negative log-likelihood is infinite at starting values")
    
    opt <- optim(start, nllh, hessian = TRUE, ..., method = method)
    if (opt$convergence > 0) 
        warning("optimization may not have succeeded")
    if(std.err) {
        tol <- .Machine$double.eps^0.5
        std.err <- qr(opt$hessian, tol = tol)
        if (std.err$rank != ncol(std.err$qr)) 
            stop("observed information matrix is singular; use std.err = FALSE")
        std.err <- diag(solve(std.err, tol = tol))
        if(any(std.err <= 0))
            stop("observed information matrix is singular; use std.err = FALSE")
        std.err <- sqrt(std.err)
    }
    else std.err <- rep(NA,l)
    names(std.err) <- nm
    list(estimate = opt$par, std.err = std.err, deviance = 2*opt$value,
         counts = opt$counts)
}

"fbvhr"<- 
function(x, start, ..., nsloc1 = NULL, nsloc2 = NULL, std.err = TRUE,
         method = "BFGS")
{
    nlbvhr <- function(loc1 = 0, scale1 = 1, shape1 = 0, loc2 = loc1,
                        scale2 = scale1, shape2 = shape1, dep)
    {
        if(any(c(scale1,scale2) < 0.01) || dep < 0.2 || dep > 10)
            return(1e6)
        if(!is.null(nsloc1)) {
            ns <- numeric(length(loc.param1))
            for(i in 1:length(ns))
                ns[i] <- get(loc.param1[i])
            loc1 <- drop(nslocmat1 %*% ns)
        }
        else loc1 <- rep(loc1, length.out = length(naind))
        if(!is.null(nsloc2)) {
            ns <- numeric(length(loc.param2))
            for(i in 1:length(ns))
                ns[i] <- get(loc.param2[i])
            loc2 <- drop(nslocmat2 %*% ns)
        }
        else loc2 <- rep(loc2, length.out = length(naind))
        if(any(naind == 2))
            m1l <- .C("nlgev",
                x.m1, n.m1,
                as.double(loc1[naind == 2]), as.double(scale1),
                as.double(shape1), dns = double(1), PACKAGE = "evd")$dns
        else m1l <- 0
        if(is.nan(m1l))
            stop("NaN returned in likelihood evaluation\n\tPlease contact me at a.stephenson@lancaster.ac.uk")
        if(any(naind == 1))
            m2l <- .C("nlgev",
                x.m2, n.m2,
                as.double(loc2[naind == 1]), as.double(scale2),
                as.double(shape2), dns = double(1), PACKAGE = "evd")$dns
        else m2l <- 0
        if(is.nan(m2l))
            stop("NaN returned in likelihood evaluation\n\tPlease contact me at a.stephenson@lancaster.ac.uk")
        bvl <- .C("nlbvhr",
                x1, x2, n,
                as.double(dep), as.double(loc1[naind == 0]), as.double(scale1),
                as.double(shape1), as.double(loc2[naind == 0]),
                as.double(scale2), as.double(shape2), dns = double(1),
                PACKAGE = "evd")$dns
        if(is.nan(bvl))
            stop("NaN returned in likelihood evaluation\n\tPlease contact me at a.stephenson@lancaster.ac.uk")
        if(any(c(m1l,m2l,bvl) == 1e6)) return(1e6)
        else return(m1l + m2l + bvl)
    }
    naind <- na.vals(x)
    if(!any(naind == 0))
        stop("`x' must have at least one complete observation")
    if(!is.null(nsloc1)) {
        if(is.vector(nsloc1))
            nsloc1 <- data.frame(trend = nsloc1)
        if(!is.data.frame(nsloc1))
            stop("`nsloc1' must be a vector or data frame")
        if(nrow(nsloc1) != nrow(x))
            stop("`nsloc1' and `x' are not compatible")
        nslocmat1 <- cbind(1,as.matrix(nsloc1))
    }
    if(!is.null(nsloc2)) {
        if(is.vector(nsloc2))
            nsloc2 <- data.frame(trend = nsloc2)
        if(!is.data.frame(nsloc2))
            stop("`nsloc2' must be a vector or data frame")
        if(nrow(nsloc2) != nrow(x))
            stop("`nsloc2' and `x' are not compatible")
        nslocmat2 <- cbind(1,as.matrix(nsloc2))
    }
    loc.param1 <- paste("loc1", c("",names(nsloc1)), sep="")
    loc.param2 <- paste("loc2", c("",names(nsloc2)), sep="")
    param <- c(loc.param1, "scale1", "shape1", loc.param2, "scale2",
               "shape2", "dep")
    if(missing(start)) {
        start <- as.list(numeric(length(param)))
        names(start) <- param
        st1 <- as.list(fgev(x[,1], std.err = FALSE, nsloc = nsloc1)$estimate)
        st2 <- as.list(fgev(x[,2], std.err = FALSE, nsloc = nsloc2)$estimate)
        start[c(loc.param1, "scale1", "shape1")] <- st1
        start[c(loc.param2, "scale2", "shape2")] <- st2
        start[["dep"]] <- 1
        start <- start[!(param %in% names(list(...)))]
    }
    if(any(!is.na(match(names(start),c("mar1","mar2"))))) {
        start <- unlist(start)
        names(start)[grep("mar1",names(start))] <- c("loc1","scale1","shape1")
        names(start)[grep("mar2",names(start))] <- c("loc2","scale2","shape2")
        start <- as.list(start)
    }
    if(!is.list(start)) 
        stop("`start' must be a named list")
    x.m1 <- as.double(x[naind == 2, 1])
    n.m1 <- as.integer(length(x.m1))
    x.m2 <- as.double(x[naind == 1, 2])
    n.m2 <- as.integer(length(x.m2))
    x <- x[naind == 0, , drop = FALSE]
    x1 <- as.double(x[,1])
    x2 <- as.double(x[,2])
    n <- as.integer(nrow(x))
    nm <- names(start)
    l <- length(nm)
    f <- c(as.list(numeric(length(loc.param1))), formals(nlbvhr)[2:3],
           as.list(numeric(length(loc.param2))), formals(nlbvhr)[5:7])
    names(f) <- param
    m <- match(nm, param)
    if(any(is.na(m))) 
        stop("`start' specifies unknown arguments")    
    formals(nlbvhr) <- c(f[m], f[-m])
    nllh <- function(p, ...) nlbvhr(p, ...)
    if(l > 1)
        body(nllh) <- parse(text = paste("nlbvhr(", paste("p[",1:l,
            "]", collapse = ", "), ", ...)"))
    fixed.param <- list(...)[names(list(...)) %in% param]
    start.arg <- c(list(p = unlist(start)), fixed.param)
    if(do.call("nllh", start.arg) == 1e6)
        warning("negative log-likelihood is infinite at starting values")
    opt <- optim(start, nllh, hessian = TRUE, ..., method = method)
    if (opt$convergence > 0) 
        warning("optimization may not have succeeded")
    if(std.err) {
        tol <- .Machine$double.eps^0.5
        std.err <- qr(opt$hessian, tol = tol)
        if (std.err$rank != ncol(std.err$qr)) 
            stop("observed information matrix is singular; use std.err = FALSE")
        std.err <- diag(solve(std.err, tol = tol))
        if(any(std.err <= 0))
            stop("observed information matrix is singular; use std.err = FALSE")
        std.err <- sqrt(std.err)
    }
    else std.err <- rep(NA,l)
    names(std.err) <- nm
    list(estimate = opt$par, std.err = std.err, deviance = 2*opt$value,
         counts = opt$counts)
}

"fbvneglog"<- 
function(x, start, ..., nsloc1 = NULL, nsloc2 = NULL, std.err = TRUE,
         method = "BFGS")
{
    nlbvneglog <- function(loc1 = 0, scale1 = 1, shape1 = 0, loc2 = loc1,
                        scale2 = scale1, shape2 = shape1, dep)
    {
        if(any(c(scale1,scale2) < 0.01) || dep < 0.05 || dep > 8)
            return(1e6)
        if(!is.null(nsloc1)) {
            ns <- numeric(length(loc.param1))
            for(i in 1:length(ns))
                ns[i] <- get(loc.param1[i])
            loc1 <- drop(nslocmat1 %*% ns)
        }
        else loc1 <- rep(loc1, length.out = length(naind))
        if(!is.null(nsloc2)) {
            ns <- numeric(length(loc.param2))
            for(i in 1:length(ns))
                ns[i] <- get(loc.param2[i])
            loc2 <- drop(nslocmat2 %*% ns)
        }
        else loc2 <- rep(loc2, length.out = length(naind))
        if(any(naind == 2))
            m1l <- .C("nlgev",
                x.m1, n.m1,
                as.double(loc1[naind == 2]), as.double(scale1),
                as.double(shape1), dns = double(1), PACKAGE = "evd")$dns
        else m1l <- 0
        if(is.nan(m1l))
            stop("NaN returned in likelihood evaluation\n\tPlease contact me at a.stephenson@lancaster.ac.uk")
        if(any(naind == 1))
            m2l <- .C("nlgev",
                x.m2, n.m2,
                as.double(loc2[naind == 1]), as.double(scale2),
                as.double(shape2), dns = double(1), PACKAGE = "evd")$dns
        else m2l <- 0
        if(is.nan(m2l))
            stop("NaN returned in likelihood evaluation\n\tPlease contact me at a.stephenson@lancaster.ac.uk")
        bvl <- .C("nlbvneglog",
                x1, x2, n,
                as.double(dep), as.double(loc1[naind == 0]), as.double(scale1),
                as.double(shape1), as.double(loc2[naind == 0]),
                as.double(scale2), as.double(shape2), dns = double(1),
                PACKAGE = "evd")$dns
        if(is.nan(bvl))
            stop("NaN returned in likelihood evaluation\n\tPlease contact me at a.stephenson@lancaster.ac.uk")
        if(any(c(m1l,m2l,bvl) == 1e6)) return(1e6)
        else return(m1l + m2l + bvl)
    }
    naind <- na.vals(x)
    if(!any(naind == 0))
        stop("`x' must have at least one complete observation")
    if(!is.null(nsloc1)) {
        if(is.vector(nsloc1))
            nsloc1 <- data.frame(trend = nsloc1)
        if(!is.data.frame(nsloc1))
            stop("`nsloc1' must be a vector or data frame")
        if(nrow(nsloc1) != nrow(x))
            stop("`nsloc1' and `x' are not compatible")
        nslocmat1 <- cbind(1,as.matrix(nsloc1))
    }
    if(!is.null(nsloc2)) {
        if(is.vector(nsloc2))
            nsloc2 <- data.frame(trend = nsloc2)
        if(!is.data.frame(nsloc2))
            stop("`nsloc2' must be a vector or data frame")
        if(nrow(nsloc2) != nrow(x))
            stop("`nsloc2' and `x' are not compatible")
        nslocmat2 <- cbind(1,as.matrix(nsloc2))
    }
    loc.param1 <- paste("loc1", c("",names(nsloc1)), sep="")
    loc.param2 <- paste("loc2", c("",names(nsloc2)), sep="")
    param <- c(loc.param1, "scale1", "shape1", loc.param2, "scale2",
               "shape2", "dep")
    if(missing(start)) {
        start <- as.list(numeric(length(param)))
        names(start) <- param
        st1 <- as.list(fgev(x[,1], std.err = FALSE, nsloc = nsloc1)$estimate)
        st2 <- as.list(fgev(x[,2], std.err = FALSE, nsloc = nsloc2)$estimate)
        start[c(loc.param1, "scale1", "shape1")] <- st1
        start[c(loc.param2, "scale2", "shape2")] <- st2
        start[["dep"]] <- 0.6
        start <- start[!(param %in% names(list(...)))]
    }
    if(any(!is.na(match(names(start),c("mar1","mar2"))))) {
        start <- unlist(start)
        names(start)[grep("mar1",names(start))] <- c("loc1","scale1","shape1")
        names(start)[grep("mar2",names(start))] <- c("loc2","scale2","shape2")
        start <- as.list(start)
    }
    if(!is.list(start)) 
        stop("`start' must be a named list")
    x.m1 <- as.double(x[naind == 2, 1])
    n.m1 <- as.integer(length(x.m1))
    x.m2 <- as.double(x[naind == 1, 2])
    n.m2 <- as.integer(length(x.m2))
    x <- x[naind == 0, , drop = FALSE]
    x1 <- as.double(x[,1])
    x2 <- as.double(x[,2])
    n <- as.integer(nrow(x))
    nm <- names(start)
    l <- length(nm)
    f <- c(as.list(numeric(length(loc.param1))), formals(nlbvneglog)[2:3],
           as.list(numeric(length(loc.param2))), formals(nlbvneglog)[5:7])
    names(f) <- param
    
    m <- match(nm, param)
    if(any(is.na(m))) 
        stop("`start' specifies unknown arguments")    
    formals(nlbvneglog) <- c(f[m], f[-m])
    nllh <- function(p, ...) nlbvneglog(p, ...)
    if(l > 1)
        body(nllh) <- parse(text = paste("nlbvneglog(", paste("p[",1:l,
            "]", collapse = ", "), ", ...)"))
    fixed.param <- list(...)[names(list(...)) %in% param]
    start.arg <- c(list(p = unlist(start)), fixed.param)
    if(do.call("nllh", start.arg) == 1e6)
        warning("negative log-likelihood is infinite at starting values")
    opt <- optim(start, nllh, hessian = TRUE, ..., method = method)
    if (opt$convergence > 0) 
        warning("optimization may not have succeeded")
    if(std.err) {
        tol <- .Machine$double.eps^0.5
        std.err <- qr(opt$hessian, tol = tol)
        if (std.err$rank != ncol(std.err$qr)) 
            stop("observed information matrix is singular; use std.err = FALSE")
        std.err <- diag(solve(std.err, tol = tol))
        if(any(std.err <= 0))
            stop("observed information matrix is singular; use std.err = FALSE")
        std.err <- sqrt(std.err)
    }
    else std.err <- rep(NA,l)
    names(std.err) <- nm
    list(estimate = opt$par, std.err = std.err, deviance = 2*opt$value,
         counts = opt$counts)
}

"fbvaneglog"<- 
function(x, start, ..., nsloc1 = NULL, nsloc2 = NULL, std.err = TRUE,
         method = "BFGS")
{
    nlbvaneglog <- function(loc1 = 0, scale1 = 1, shape1 = 0, loc2 = loc1,
                scale2 = scale1, shape2 = shape1, asy1 = 1, asy2 = 1, dep)
    {
        if(any(c(scale1,scale2) < 0.01) || any(c(asy1,asy2) > 1) ||
           any(c(asy1,asy2) < 0.001) || dep < 0.05 || dep > 8)
            return(1e6)
        if(!is.null(nsloc1)) {
            ns <- numeric(length(loc.param1))
            for(i in 1:length(ns))
                ns[i] <- get(loc.param1[i])
            loc1 <- drop(nslocmat1 %*% ns)
        }
        else loc1 <- rep(loc1, length.out = length(naind))
        if(!is.null(nsloc2)) {
            ns <- numeric(length(loc.param2))
            for(i in 1:length(ns))
                ns[i] <- get(loc.param2[i])
            loc2 <- drop(nslocmat2 %*% ns)
        }
        else loc2 <- rep(loc2, length.out = length(naind))
        if(any(naind == 2))
            m1l <- .C("nlgev",
                x.m1, n.m1,
                as.double(loc1[naind == 2]), as.double(scale1),
                as.double(shape1), dns = double(1), PACKAGE = "evd")$dns
        else m1l <- 0
        if(is.nan(m1l))
            stop("NaN returned in likelihood evaluation\n\tPlease contact me at a.stephenson@lancaster.ac.uk")
        if(any(naind == 1))
            m2l <- .C("nlgev",
                x.m2, n.m2,
                as.double(loc2[naind == 1]), as.double(scale2),
                as.double(shape2), dns = double(1), PACKAGE = "evd")$dns
        else m2l <- 0
        if(is.nan(m2l))
            stop("NaN returned in likelihood evaluation\n\tPlease contact me at a.stephenson@lancaster.ac.uk")
        bvl <- .C("nlbvaneglog",
                x1, x2, n,
                as.double(dep), as.double(asy1), as.double(asy2),
                as.double(loc1[naind == 0]), as.double(scale1),
                as.double(shape1), as.double(loc2[naind == 0]),
                as.double(scale2), as.double(shape2), dns = double(1),
                PACKAGE = "evd")$dns
        if(is.nan(bvl))
            stop("NaN returned in likelihood evaluation\n\tPlease contact me at a.stephenson@lancaster.ac.uk")
        
        if(any(c(m1l,m2l,bvl) == 1e6)) return(1e6)
        else return(m1l + m2l + bvl)
    }
    naind <- na.vals(x)
    if(!any(naind == 0))
        stop("`x' must have at least one complete observation")
    if(!is.null(nsloc1)) {
        if(is.vector(nsloc1))
            nsloc1 <- data.frame(trend = nsloc1)
        if(!is.data.frame(nsloc1))
            stop("`nsloc1' must be a vector or data frame")
        if(nrow(nsloc1) != nrow(x))
            stop("`nsloc1' and `x' are not compatible")
        nslocmat1 <- cbind(1,as.matrix(nsloc1))
    }
    if(!is.null(nsloc2)) {
        if(is.vector(nsloc2))
            nsloc2 <- data.frame(trend = nsloc2)
        if(!is.data.frame(nsloc2))
            stop("`nsloc2' must be a vector or data frame")
        if(nrow(nsloc2) != nrow(x))
            stop("`nsloc2' and `x' are not compatible")
        nslocmat2 <- cbind(1,as.matrix(nsloc2))
    }
    loc.param1 <- paste("loc1", c("",names(nsloc1)), sep="")
    loc.param2 <- paste("loc2", c("",names(nsloc2)), sep="")
    param <- c(loc.param1, "scale1", "shape1", loc.param2, "scale2",
               "shape2", "asy1", "asy2", "dep")
    if(missing(start)) {
        start <- as.list(numeric(length(param)))
        names(start) <- param
        st1 <- as.list(fgev(x[,1], std.err = FALSE, nsloc = nsloc1)$estimate)
        st2 <- as.list(fgev(x[,2], std.err = FALSE, nsloc = nsloc2)$estimate)
        start[c(loc.param1, "scale1", "shape1")] <- st1
        start[c(loc.param2, "scale2", "shape2")] <- st2
        start[["asy1"]] <- 0.75
        start[["asy2"]] <- 0.75
        start[["dep"]] <- 0.8
        start <- start[!(param %in% names(list(...)))]
    }
    if(any(!is.na(match(names(start),c("asy","mar1","mar2"))))) {
        start <- unlist(start)
        names(start)[grep("mar1",names(start))] <- c("loc1","scale1","shape1")
        names(start)[grep("mar2",names(start))] <- c("loc2","scale2","shape2")
        start <- as.list(start)
    }
    if(!is.list(start)) 
        stop("`start' must be a named list")
    x.m1 <- as.double(x[naind == 2, 1])
    n.m1 <- as.integer(length(x.m1))
    x.m2 <- as.double(x[naind == 1, 2])
    n.m2 <- as.integer(length(x.m2))
    x <- x[naind == 0, , drop = FALSE]
    x1 <- as.double(x[,1])
    x2 <- as.double(x[,2])
    n <- as.integer(nrow(x))
    nm <- names(start)
    l <- length(nm)
    f <- c(as.list(numeric(length(loc.param1))), formals(nlbvaneglog)[2:3],
           as.list(numeric(length(loc.param2))), formals(nlbvaneglog)[5:9])
    names(f) <- param
    m <- match(nm, param)
    if(any(is.na(m))) 
        stop("`start' specifies unknown arguments")    
    formals(nlbvaneglog) <- c(f[m], f[-m])
    nllh <- function(p, ...) nlbvaneglog(p, ...)
    if(l > 1)
        body(nllh) <- parse(text = paste("nlbvaneglog(", paste("p[",1:l,
            "]", collapse = ", "), ", ...)"))
    fixed.param <- list(...)[names(list(...)) %in% param]
    start.arg <- c(list(p = unlist(start)), fixed.param)
    if(do.call("nllh", start.arg) == 1e6)
        warning("negative log-likelihood is infinite at starting values")
    
    opt <- optim(start, nllh, hessian = TRUE, ..., method = method)
    if (opt$convergence > 0) 
        warning("optimization may not have succeeded")
    if(std.err) {
        tol <- .Machine$double.eps^0.5
        std.err <- qr(opt$hessian, tol = tol)
        if (std.err$rank != ncol(std.err$qr)) 
            stop("observed information matrix is singular; use std.err = FALSE")
        std.err <- diag(solve(std.err, tol = tol))
        if(any(std.err <= 0))
            stop("observed information matrix is singular; use std.err = FALSE")
        std.err <- sqrt(std.err)
    }
    else std.err <- rep(NA,l)
    names(std.err) <- nm
    list(estimate = opt$par, std.err = std.err, deviance = 2*opt$value,
         counts = opt$counts)
}

"fbvbilog"<- 
function(x, start, ..., nsloc1 = NULL, nsloc2 = NULL, std.err = TRUE,
         method = "BFGS")
{
    nlbvbilog <- function(loc1 = 0, scale1 = 1, shape1 = 0, loc2 = loc1,
                        scale2 = scale1, shape2 = shape1, alpha, beta)
    {
        if(any(c(scale1,scale2) < 0.01) || any(c(alpha,beta) < 0.1) ||
           any(c(alpha,beta) > 0.999))
            return(1e6)
        if(!is.null(nsloc1)) {
            ns <- numeric(length(loc.param1))
            for(i in 1:length(ns))
                ns[i] <- get(loc.param1[i])
            loc1 <- drop(nslocmat1 %*% ns)
        }
        else loc1 <- rep(loc1, length.out = length(naind))
        if(!is.null(nsloc2)) {
            ns <- numeric(length(loc.param2))
            for(i in 1:length(ns))
                ns[i] <- get(loc.param2[i])
            loc2 <- drop(nslocmat2 %*% ns)
        }
        else loc2 <- rep(loc2, length.out = length(naind))
        if(any(naind == 2))
            m1l <- .C("nlgev",
                x.m1, n.m1,
                as.double(loc1[naind == 2]), as.double(scale1),
                as.double(shape1), dns = double(1), PACKAGE = "evd")$dns
        else m1l <- 0
        if(is.nan(m1l))
            stop("NaN returned in likelihood evaluation\n\tPlease contact me at a.stephenson@lancaster.ac.uk")
        if(any(naind == 1))
            m2l <- .C("nlgev",
                x.m2, n.m2,
                as.double(loc2[naind == 1]), as.double(scale2),
                as.double(shape2), dns = double(1), PACKAGE = "evd")$dns
        else m2l <- 0
        if(is.nan(m2l))
            stop("NaN returned in likelihood evaluation\n\tPlease contact me at a.stephenson@lancaster.ac.uk")
        bvl <- .C("nlbvbilog",
                x1, x2, n,
                as.double(alpha), as.double(beta), as.double(loc1[naind == 0]),
                as.double(scale1), as.double(shape1),
                as.double(loc2[naind == 0]), as.double(scale2),
                as.double(shape2), dns = double(1), PACKAGE = "evd")$dns
        if(is.nan(bvl))
            stop("NaN returned in likelihood evaluation\n\tPlease contact me at a.stephenson@lancaster.ac.uk")
        if(any(c(m1l,m2l,bvl) >= 1e6)) return(1e6)
        else return(m1l + m2l + bvl)
    }
    naind <- na.vals(x)
    if(!any(naind == 0))
        stop("`x' must have at least one complete observation")
    if(!is.null(nsloc1)) {
        if(is.vector(nsloc1))
            nsloc1 <- data.frame(trend = nsloc1)
        if(!is.data.frame(nsloc1))
            stop("`nsloc1' must be a vector or data frame")
        if(nrow(nsloc1) != nrow(x))
            stop("`nsloc1' and `x' are not compatible")
        nslocmat1 <- cbind(1,as.matrix(nsloc1))
    }
    if(!is.null(nsloc2)) {
        if(is.vector(nsloc2))
            nsloc2 <- data.frame(trend = nsloc2)
        if(!is.data.frame(nsloc2))
            stop("`nsloc2' must be a vector or data frame")
        if(nrow(nsloc2) != nrow(x))
            stop("`nsloc2' and `x' are not compatible")
        nslocmat2 <- cbind(1,as.matrix(nsloc2))
    }
    loc.param1 <- paste("loc1", c("",names(nsloc1)), sep="")
    loc.param2 <- paste("loc2", c("",names(nsloc2)), sep="")
    param <- c(loc.param1, "scale1", "shape1", loc.param2, "scale2",
               "shape2", "alpha", "beta")
    if(missing(start)) {
        start <- as.list(numeric(length(param)))
        names(start) <- param
        st1 <- as.list(fgev(x[,1], std.err = FALSE, nsloc = nsloc1)$estimate)
        st2 <- as.list(fgev(x[,2], std.err = FALSE, nsloc = nsloc2)$estimate)
        start[c(loc.param1, "scale1", "shape1")] <- st1
        start[c(loc.param2, "scale2", "shape2")] <- st2
        start[["alpha"]] <- 0.75
        start[["beta"]] <- 0.75
        start <- start[!(param %in% names(list(...)))]
    }
    if(any(!is.na(match(names(start),c("mar1","mar2"))))) {
        start <- unlist(start)
        names(start)[grep("mar1",names(start))] <- c("loc1","scale1","shape1")
        names(start)[grep("mar2",names(start))] <- c("loc2","scale2","shape2")
        start <- as.list(start)
    }
    if(!is.list(start)) 
        stop("`start' must be a named list")
    x.m1 <- as.double(x[naind == 2, 1])
    n.m1 <- as.integer(length(x.m1))
    x.m2 <- as.double(x[naind == 1, 2])
    n.m2 <- as.integer(length(x.m2))
    x <- x[naind == 0, , drop = FALSE]
    x1 <- as.double(x[,1])
    x2 <- as.double(x[,2])
    n <- as.integer(nrow(x))
    nm <- names(start)
    l <- length(nm)
    f <- c(as.list(numeric(length(loc.param1))), formals(nlbvbilog)[2:3],
           as.list(numeric(length(loc.param2))), formals(nlbvbilog)[5:8])
    names(f) <- param
    m <- match(nm, param)
    if(any(is.na(m))) 
        stop("`start' specifies unknown arguments")    
    formals(nlbvbilog) <- c(f[m], f[-m])
    nllh <- function(p, ...) nlbvbilog(p, ...)
    if(l > 1)
        body(nllh) <- parse(text = paste("nlbvbilog(", paste("p[",1:l,
            "]", collapse = ", "), ", ...)"))
    fixed.param <- list(...)[names(list(...)) %in% param]
    start.arg <- c(list(p = unlist(start)), fixed.param)
    if(do.call("nllh", start.arg) == 1e6)
        warning("negative log-likelihood is infinite at starting values")
    
    opt <- optim(start, nllh, hessian = TRUE, ..., method = method)
    if (opt$convergence > 0) 
        warning("optimization may not have succeeded")
    if(std.err) {
        tol <- .Machine$double.eps^0.5
        std.err <- qr(opt$hessian, tol = tol)
        if (std.err$rank != ncol(std.err$qr)) 
            stop("observed information matrix is singular; use std.err = FALSE")
        std.err <- diag(solve(std.err, tol = tol))
        if(any(std.err <= 0))
            stop("observed information matrix is singular; use std.err = FALSE")
        std.err <- sqrt(std.err)
    }
    else std.err <- rep(NA,l)
    names(std.err) <- nm
    list(estimate = opt$par, std.err = std.err, deviance = 2*opt$value,
         counts = opt$counts)
}

"fbvnegbilog"<- 
function(x, start, ..., nsloc1 = NULL, nsloc2 = NULL, std.err = TRUE,
         method = "BFGS")
{
    nlbvnegbilog <- function(loc1 = 0, scale1 = 1, shape1 = 0, loc2 = loc1,
                        scale2 = scale1, shape2 = shape1, alpha, beta)
    {
        if(any(c(scale1,scale2) < 0.01) || any(c(alpha,beta) < 0.1) ||
           any(c(alpha,beta) > 20))
            return(1e6)
        if(!is.null(nsloc1)) {
            ns <- numeric(length(loc.param1))
            for(i in 1:length(ns))
                ns[i] <- get(loc.param1[i])
            loc1 <- drop(nslocmat1 %*% ns)
        }
        else loc1 <- rep(loc1, length.out = length(naind))
        if(!is.null(nsloc2)) {
            ns <- numeric(length(loc.param2))
            for(i in 1:length(ns))
                ns[i] <- get(loc.param2[i])
            loc2 <- drop(nslocmat2 %*% ns)
        }
        else loc2 <- rep(loc2, length.out = length(naind))
        if(any(naind == 2))
            m1l <- .C("nlgev",
                x.m1, n.m1,
                as.double(loc1[naind == 2]), as.double(scale1),
                as.double(shape1), dns = double(1), PACKAGE = "evd")$dns
        else m1l <- 0
        if(is.nan(m1l))
            stop("NaN returned in likelihood evaluation\n\tPlease contact me at a.stephenson@lancaster.ac.uk")
        if(any(naind == 1))
            m2l <- .C("nlgev",
                x.m2, n.m2,
                as.double(loc2[naind == 1]), as.double(scale2),
                as.double(shape2), dns = double(1), PACKAGE = "evd")$dns
        else m2l <- 0
        if(is.nan(m2l))
            stop("NaN returned in likelihood evaluation\n\tPlease contact me at a.stephenson@lancaster.ac.uk")
        bvl <- .C("nlbvnegbilog",
                x1, x2, n,
                as.double(alpha), as.double(beta), as.double(loc1[naind == 0]),
                as.double(scale1), as.double(shape1),
                as.double(loc2[naind == 0]), as.double(scale2),
                as.double(shape2), dns = double(1), PACKAGE = "evd")$dns
        if(is.nan(bvl))
            stop("NaN returned in likelihood evaluation\n\tPlease contact me at a.stephenson@lancaster.ac.uk")
        if(any(c(m1l,m2l,bvl) >= 1e6)) return(1e6)
        else return(m1l + m2l + bvl)
    }
    naind <- na.vals(x)
    if(!any(naind == 0))
        stop("`x' must have at least one complete observation")
    if(!is.null(nsloc1)) {
        if(is.vector(nsloc1))
            nsloc1 <- data.frame(trend = nsloc1)
        if(!is.data.frame(nsloc1))
            stop("`nsloc1' must be a vector or data frame")
        if(nrow(nsloc1) != nrow(x))
            stop("`nsloc1' and `x' are not compatible")
        nslocmat1 <- cbind(1,as.matrix(nsloc1))
    }
    if(!is.null(nsloc2)) {
        if(is.vector(nsloc2))
            nsloc2 <- data.frame(trend = nsloc2)
        if(!is.data.frame(nsloc2))
            stop("`nsloc2' must be a vector or data frame")
        if(nrow(nsloc2) != nrow(x))
            stop("`nsloc2' and `x' are not compatible")
        nslocmat2 <- cbind(1,as.matrix(nsloc2))
    }
    loc.param1 <- paste("loc1", c("",names(nsloc1)), sep="")
    loc.param2 <- paste("loc2", c("",names(nsloc2)), sep="")
    param <- c(loc.param1, "scale1", "shape1", loc.param2, "scale2",
               "shape2", "alpha", "beta")
    if(missing(start)) {
        start <- as.list(numeric(length(param)))
        names(start) <- param
        st1 <- as.list(fgev(x[,1], std.err = FALSE, nsloc = nsloc1)$estimate)
        st2 <- as.list(fgev(x[,2], std.err = FALSE, nsloc = nsloc2)$estimate)
        start[c(loc.param1, "scale1", "shape1")] <- st1
        start[c(loc.param2, "scale2", "shape2")] <- st2
        start[["alpha"]] <- 1/0.6
        start[["beta"]] <- 1/0.6
        start <- start[!(param %in% names(list(...)))]
    }
    if(any(!is.na(match(names(start),c("mar1","mar2"))))) {
        start <- unlist(start)
        names(start)[grep("mar1",names(start))] <- c("loc1","scale1","shape1")
        names(start)[grep("mar2",names(start))] <- c("loc2","scale2","shape2")
        start <- as.list(start)
    }
    if(!is.list(start)) 
        stop("`start' must be a named list")
    x.m1 <- as.double(x[naind == 2, 1])
    n.m1 <- as.integer(length(x.m1))
    x.m2 <- as.double(x[naind == 1, 2])
    n.m2 <- as.integer(length(x.m2))
    x <- x[naind == 0, , drop = FALSE]
    x1 <- as.double(x[,1])
    x2 <- as.double(x[,2])
    n <- as.integer(nrow(x))
    nm <- names(start)
    l <- length(nm)
    f <- c(as.list(numeric(length(loc.param1))), formals(nlbvnegbilog)[2:3],
           as.list(numeric(length(loc.param2))), formals(nlbvnegbilog)[5:8])
    names(f) <- param
    m <- match(nm, param)
    if(any(is.na(m))) 
        stop("`start' specifies unknown arguments")    
    formals(nlbvnegbilog) <- c(f[m], f[-m])
    nllh <- function(p, ...) nlbvnegbilog(p, ...)
    if(l > 1)
        body(nllh) <- parse(text = paste("nlbvnegbilog(", paste("p[",1:l,
            "]", collapse = ", "), ", ...)"))
    fixed.param <- list(...)[names(list(...)) %in% param]
    start.arg <- c(list(p = unlist(start)), fixed.param)
    if(do.call("nllh", start.arg) == 1e6)
        warning("negative log-likelihood is infinite at starting values")
    
    opt <- optim(start, nllh, hessian = TRUE, ..., method = method)
    if (opt$convergence > 0) 
        warning("optimization may not have succeeded")
    if(std.err) {
        tol <- .Machine$double.eps^0.5
        std.err <- qr(opt$hessian, tol = tol)
        if (std.err$rank != ncol(std.err$qr)) 
            stop("observed information matrix is singular; use std.err = FALSE")
        std.err <- diag(solve(std.err, tol = tol))
        if(any(std.err <= 0))
            stop("observed information matrix is singular; use std.err = FALSE")
        std.err <- sqrt(std.err)
    }
    else std.err <- rep(NA,l)
    names(std.err) <- nm
    list(estimate = opt$par, std.err = std.err, deviance = 2*opt$value,
         counts = opt$counts)
}

"fbvct"<- 
function(x, start, ..., nsloc1 = NULL, nsloc2 = NULL, std.err = TRUE,
         method = "BFGS")
{
    nlbvct <- function(loc1 = 0, scale1 = 1, shape1 = 0, loc2 = loc1,
                        scale2 = scale1, shape2 = shape1, alpha, beta)
    {
        if(any(c(scale1,scale2) < 0.01) || any(c(alpha,beta) < 0.001) ||
           any(c(alpha,beta) > 30))
            return(1e6)
        if(!is.null(nsloc1)) {
            ns <- numeric(length(loc.param1))
            for(i in 1:length(ns))
                ns[i] <- get(loc.param1[i])
            loc1 <- drop(nslocmat1 %*% ns)
        }
        else loc1 <- rep(loc1, length.out = length(naind))
        if(!is.null(nsloc2)) {
            ns <- numeric(length(loc.param2))
            for(i in 1:length(ns))
                ns[i] <- get(loc.param2[i])
            loc2 <- drop(nslocmat2 %*% ns)
        }
        else loc2 <- rep(loc2, length.out = length(naind))
        if(any(naind == 2))
            m1l <- .C("nlgev",
                x.m1, n.m1,
                as.double(loc1[naind == 2]), as.double(scale1),
                as.double(shape1), dns = double(1), PACKAGE = "evd")$dns
        else m1l <- 0
        if(is.nan(m1l))
            stop("NaN returned in likelihood evaluation\n\tPlease contact me at a.stephenson@lancaster.ac.uk")
        if(any(naind == 1))
            m2l <- .C("nlgev",
                x.m2, n.m2,
                as.double(loc2[naind == 1]), as.double(scale2),
                as.double(shape2), dns = double(1), PACKAGE = "evd")$dns
        else m2l <- 0
        if(is.nan(m2l))
            stop("NaN returned in likelihood evaluation\n\tPlease contact me at a.stephenson@lancaster.ac.uk")
        bvl <- .C("nlbvct",
                x1, x2, n,
                as.double(alpha), as.double(beta), as.double(loc1[naind == 0]),
                as.double(scale1), as.double(shape1),
                as.double(loc2[naind == 0]), as.double(scale2),
                as.double(shape2), dns = double(1), PACKAGE = "evd")$dns
        if(is.nan(bvl))
            stop("NaN returned in likelihood evaluation\n\tPlease contact me at a.stephenson@lancaster.ac.uk")
        if(any(c(m1l,m2l,bvl) >= 1e6)) return(1e6)
        else return(m1l + m2l + bvl)
    }
    naind <- na.vals(x)
    if(!any(naind == 0))
        stop("`x' must have at least one complete observation")
    if(!is.null(nsloc1)) {
        if(is.vector(nsloc1))
            nsloc1 <- data.frame(trend = nsloc1)
        if(!is.data.frame(nsloc1))
            stop("`nsloc1' must be a vector or data frame")
        if(nrow(nsloc1) != nrow(x))
            stop("`nsloc1' and `x' are not compatible")
        nslocmat1 <- cbind(1,as.matrix(nsloc1))
    }
    if(!is.null(nsloc2)) {
        if(is.vector(nsloc2))
            nsloc2 <- data.frame(trend = nsloc2)
        if(!is.data.frame(nsloc2))
            stop("`nsloc2' must be a vector or data frame")
        if(nrow(nsloc2) != nrow(x))
            stop("`nsloc2' and `x' are not compatible")
        nslocmat2 <- cbind(1,as.matrix(nsloc2))
    }
    loc.param1 <- paste("loc1", c("",names(nsloc1)), sep="")
    loc.param2 <- paste("loc2", c("",names(nsloc2)), sep="")
    param <- c(loc.param1, "scale1", "shape1", loc.param2, "scale2",
               "shape2", "alpha", "beta")
    if(missing(start)) {
        start <- as.list(numeric(length(param)))
        names(start) <- param
        st1 <- as.list(fgev(x[,1], std.err = FALSE, nsloc = nsloc1)$estimate)
        st2 <- as.list(fgev(x[,2], std.err = FALSE, nsloc = nsloc2)$estimate)
        start[c(loc.param1, "scale1", "shape1")] <- st1
        start[c(loc.param2, "scale2", "shape2")] <- st2
        start[["alpha"]] <- 0.6
        start[["beta"]] <- 0.6
        start <- start[!(param %in% names(list(...)))]
    }
    if(any(!is.na(match(names(start),c("mar1","mar2"))))) {
        start <- unlist(start)
        names(start)[grep("mar1",names(start))] <- c("loc1","scale1","shape1")
        names(start)[grep("mar2",names(start))] <- c("loc2","scale2","shape2")
        start <- as.list(start)
    }
    if(!is.list(start)) 
        stop("`start' must be a named list")
    x.m1 <- as.double(x[naind == 2, 1])
    n.m1 <- as.integer(length(x.m1))
    x.m2 <- as.double(x[naind == 1, 2])
    n.m2 <- as.integer(length(x.m2))
    x <- x[naind == 0, , drop = FALSE]
    x1 <- as.double(x[,1])
    x2 <- as.double(x[,2])
    n <- as.integer(nrow(x))
    nm <- names(start)
    l <- length(nm)
    f <- c(as.list(numeric(length(loc.param1))), formals(nlbvct)[2:3],
           as.list(numeric(length(loc.param2))), formals(nlbvct)[5:8])
    names(f) <- param
    m <- match(nm, param)
    if(any(is.na(m))) 
        stop("`start' specifies unknown arguments")    
    formals(nlbvct) <- c(f[m], f[-m])
    nllh <- function(p, ...) nlbvct(p, ...)
    if(l > 1)
        body(nllh) <- parse(text = paste("nlbvct(", paste("p[",1:l,
            "]", collapse = ", "), ", ...)"))
    fixed.param <- list(...)[names(list(...)) %in% param]
    start.arg <- c(list(p = unlist(start)), fixed.param)
    if(do.call("nllh", start.arg) == 1e6)
        warning("negative log-likelihood is infinite at starting values")
    
    opt <- optim(start, nllh, hessian = TRUE, ..., method = method)
    if (opt$convergence > 0) 
        warning("optimization may not have succeeded")
    if(std.err) {
        tol <- .Machine$double.eps^0.5
        std.err <- qr(opt$hessian, tol = tol)
        if (std.err$rank != ncol(std.err$qr)) 
            stop("observed information matrix is singular; use std.err = FALSE")
        std.err <- diag(solve(std.err, tol = tol))
        if(any(std.err <= 0))
            stop("observed information matrix is singular; use std.err = FALSE")
        std.err <- sqrt(std.err)
    }
    else std.err <- rep(NA,l)
    names(std.err) <- nm
    list(estimate = opt$par, std.err = std.err, deviance = 2*opt$value,
         counts = opt$counts)
}

"fbvall"<-
function(x, nsloc1 = NULL, nsloc2 = NULL, omit = NULL, boxcon = TRUE, std.err = TRUE, orderby = c("AIC","BIC","SC"), control = list(maxit = 250))
{
    if(!is.null(nsloc1)) {
        if(is.vector(nsloc1))
            nsloc1 <- data.frame(trend = nsloc1)
        if(!is.data.frame(nsloc1))
            stop("`nsloc1' must be a vector or data frame")
        if(nrow(nsloc1) != nrow(x))
            stop("`nsloc1' and `x' are not compatible")
    }
    if(!is.null(nsloc2)) {
        if(is.vector(nsloc2))
            nsloc2 <- data.frame(trend = nsloc2)
        if(!is.data.frame(nsloc2))
            stop("`nsloc2' must be a vector or data frame")
        if(nrow(nsloc2) != nrow(x))
            stop("`nsloc2' and `x' are not compatible")
    }
    loc.param1 <- paste("loc1", c("",names(nsloc1)), sep="")
    loc.param2 <- paste("loc2", c("",names(nsloc2)), sep="")
    mparam <- c(loc.param1, "scale1", "shape1", loc.param2, "scale2",
               "shape2")
    nmp <- length(mparam)
    param <- c(mparam, "asy1/alpha", "asy2/beta", "dep")
    models <- c("log","alog","hr","neglog","aneglog","bilog","negbilog","ct")
    include <- !(models %in% omit)
    models <- models[include]
    std.err <- rep(std.err, length.out = length(models))
    npars <- (nmp + c(1,3,1,1,3,2,2,2))[include]
    lw <- rep(-Inf, nmp)
    up <- rep(Inf, nmp)
    n <- sum(na.vals(x) == 0)
    estimate <- std.errors <- matrix(NA, nrow=nmp+3, ncol = length(models),
                                     dimnames=list(param,models))
    criteria <- matrix(NA, nrow = 3, ncol = length(models),
                       dimnames=list(c("AIC","SC","BIC"),models))
    dm <- matrix(0, nrow = 3, ncol = length(models),
                 dimnames=list(c("dep","intdep","intasy"),models))
    if(boxcon) method <- "L-BFGS-B"
    else method <- "BFGS"
    lower <- -Inf
    upper <- Inf
    if(include[1]) {
    if(boxcon) {
        lower <- c(lw, 0.05)            
        upper <- c(up, 1)
    }
    fbvlog <- fbvlog(x = x, nsloc1 = nsloc1, nsloc2 = nsloc2,
              std.err = std.err[1], method = method,
              lower = lower, upper = upper, control = control)
    estimate[c(1:nmp,nmp+3),1] <- fbvlog$estimate
    std.errors[c(1:nmp,nmp+3),1] <- fbvlog$std.err
    dm[1,1] <- 2*(1-abvlog(dep = fbvlog$estimate["dep"]))
    dm[2,1] <- 4 * integrate(function(x)
        1-abvlog(x, dep = fbvlog$estimate["dep"]), 0, 1)$value
    print("log")
    }
    if(include[2]) {
    if(boxcon) {
        lower <- c(lw, 0, 0, 0.05)            
        upper <- c(up, 1, 1, 1)
    }
    fbvalog <- fbvalog(x = x, nsloc1 = nsloc1, nsloc2 = nsloc2,
               std.err = std.err[sum(include[1:2])], method = method,
               lower = lower, upper = upper, control = control)
    estimate[,sum(include[1:2])] <- fbvalog$estimate
    std.errors[,sum(include[1:2])] <- fbvalog$std.err
    asym <- c(fbvalog$estimate["asy1"],fbvalog$estimate["asy2"])
    depnd <- fbvalog$estimate["dep"]
    dm[1,sum(include[1:2])] <- 2*(1-abvalog(dep = depnd, asy = asym))
    depfn <- function(x) 
        abvalog(x, dep = depnd, asy = asym) -
          abvalog(x, dep = depnd, asy = c(asym[2],asym[1]))
    dm[2,sum(include[1:2])] <- 4 * integrate(function(x)
        1-abvalog(x, dep = depnd, asy = asym), 0, 1)$value
    dm[3,sum(include[1:2])] <- 4*integrate(depfn,0,.5)$value / (3-2*sqrt(2))
    print("alog")
    }
    if(include[3]) {
    if(boxcon) {
        lower <- c(lw, 0.2)            
        upper <- c(up, 10)
    }
    fbvhr <- fbvhr(x = x, nsloc1 = nsloc1, nsloc2 = nsloc2,
             std.err = std.err[sum(include[1:3])], method = method,
             lower = lower, upper = upper, control = control)
    estimate[c(1:nmp,nmp+3),sum(include[1:3])] <- fbvhr$estimate
    std.errors[c(1:nmp,nmp+3),sum(include[1:3])] <- fbvhr$std.err
    dm[1,sum(include[1:3])] <- 2*(1-abvhr(dep = fbvhr$estimate["dep"]))
    dm[2,sum(include[1:3])] <- 4 * integrate(function(x)
        1-abvhr(x, dep = fbvhr$estimate["dep"]), 0, 1)$value
    print("hr")
    }
    if(include[4]) {
    if(boxcon) {
        lower <- c(lw, 0.05)            
        upper <- c(up, 8)
    }
    fbvneglog <- fbvneglog(x = x, nsloc1 = nsloc1, nsloc2 = nsloc2,
                 std.err = std.err[sum(include[1:4])], method = method,
                 lower = lower, upper = upper, control = control)
    estimate[c(1:nmp,nmp+3),sum(include[1:4])] <- fbvneglog$estimate
    std.errors[c(1:nmp,nmp+3),sum(include[1:4])] <- fbvneglog$std.err
    dm[1,sum(include[1:4])] <- 2*(1-abvneglog(dep = fbvneglog$estimate["dep"]))
    dm[2,sum(include[1:4])] <- 4 * integrate(function(x)
        1-abvneglog(x, dep = fbvneglog$estimate["dep"]), 0, 1)$value
    print("neglog")
    }
    if(include[5]) {
    if(boxcon) {
        lower <- c(lw, 0.001, 0.001, 0.05)            
        upper <- c(up, 1, 1, 8)
    }
    fbvaneglog <- fbvaneglog(x = x, nsloc1 = nsloc1, nsloc2 = nsloc2,
                  std.err = std.err[sum(include[1:5])], method = method,
                  lower = lower, upper = upper, control = control)
    estimate[,sum(include[1:5])] <- fbvaneglog$estimate
    std.errors[,sum(include[1:5])] <- fbvaneglog$std.err
    asym <- c(fbvaneglog$estimate["asy1"],fbvaneglog$estimate["asy2"])
    depnd <- fbvaneglog$estimate["dep"]
    dm[1,sum(include[1:5])] <- 2*(1-abvaneglog(dep = depnd, asy = asym))
    depfn <- function(x) 
        abvaneglog(x, dep = depnd, asy = asym) -
          abvaneglog(x, dep = depnd, asy = c(asym[2],asym[1])) 
    dm[2,sum(include[1:5])] <- 4 * integrate(function(x)
        1-abvaneglog(x, dep = depnd, asy = asym), 0, 1)$value
    dm[3,sum(include[1:5])] <- 4*integrate(depfn,0,.5)$value / (3-2*sqrt(2))
    print("aneglog")
    }
    if(include[6]) {
    if(boxcon) {
        lower <- c(lw, 0.1, 0.1)            
        upper <- c(up, 0.999, 0.999)
    }
    fbvbilog <- fbvbilog(x = x, nsloc1 = nsloc1, nsloc2 = nsloc2,
               std.err = std.err[sum(include[1:6])], method = method,
               lower = lower, upper = upper, control = control)
    estimate[1:(nmp+2),sum(include[1:6])] <- fbvbilog$estimate
    std.errors[1:(nmp+2),sum(include[1:6])] <- fbvbilog$std.err
    alpha <- fbvbilog$estimate["alpha"]
    beta <- fbvbilog$estimate["beta"]
    dm[1,sum(include[1:6])] <- 2*(1-abvbilog(alpha = alpha, beta = beta))
    depfn <- function(x) 
        abvbilog(x, alpha = alpha, beta = beta) -
          abvbilog(x, alpha = beta, beta = alpha) 
    dm[2,sum(include[1:6])] <- 4 * integrate(function(x)
        1-abvbilog(x, alpha = alpha, beta = beta), 0, 1)$value
    dm[3,sum(include[1:6])] <- 4*integrate(depfn,0,.5)$value / (3-2*sqrt(2))
    print("bilog")
    }
    if(include[7]) {
    if(boxcon) {
        lower <- c(lw, 0.1, 0.1)             
        upper <- c(up, 20, 20) 
    }
    fbvnegbilog <- fbvnegbilog(x = x, nsloc1 = nsloc1, nsloc2 = nsloc2,
                  std.err = std.err[sum(include[1:7])], method = method,
                  lower = lower, upper = upper, control = control)
    estimate[1:(nmp+2),sum(include[1:7])] <- fbvnegbilog$estimate
    std.errors[1:(nmp+2),sum(include[1:7])] <- fbvnegbilog$std.err
    alpha <- fbvnegbilog$estimate["alpha"]
    beta <- fbvnegbilog$estimate["beta"]
    dm[1,sum(include[1:7])] <- 2*(1-abvnegbilog(alpha = alpha, beta = beta))
    depfn <- function(x) 
        abvnegbilog(x, alpha = alpha, beta = beta) -
          abvnegbilog(x, alpha = beta, beta = alpha) 
    dm[2,sum(include[1:7])] <- 4 * integrate(function(x)
        1-abvnegbilog(x, alpha = alpha, beta = beta), 0, 1)$value
    dm[3,sum(include[1:7])] <- 4*integrate(depfn,0,.5)$value / (3-2*sqrt(2))
    print("negbilog")
    }
    if(include[8]) {
    if(boxcon) {
        lower <- c(lw, 0.001, 0.001)           
        upper <- c(up, 30)
    }
    fbvct <- fbvct(x = x, nsloc1 = nsloc1, nsloc2 = nsloc2,
                   std.err = std.err[sum(include[1:8])], method = method,
                   lower = lower, upper = upper, control = control)
    estimate[1:(nmp+2),sum(include[1:8])] <- fbvct$estimate
    std.errors[1:(nmp+2),sum(include[1:8])] <- fbvct$std.err
    alpha <- fbvct$estimate["alpha"]
    beta <- fbvct$estimate["beta"]
    dm[1,sum(include[1:8])] <- 2*(1-abvct(alpha = alpha, beta = beta))
    depfn <- function(x) 
        abvct(x, alpha = alpha, beta = beta) -
          abvct(x, alpha = beta, beta = alpha) 
    dm[2,sum(include[1:8])] <- 4 * integrate(function(x)
        1-abvct(x, alpha = alpha, beta = beta), 0, 1)$value
    dm[3,sum(include[1:8])] <- 4*integrate(depfn,0,.5)$value / (3-2*sqrt(2))
    print("ct")
    }
    devs <- numeric(length(models))
    fits <- paste("fbv",models,sep="")
    for(i in 1:length(models)) 
        devs[i] <- get(fits[i])$deviance 
    criteria[1,] <- devs + 2 * npars
    criteria[2,] <- devs + log(n) * npars
    criteria[3,] <- devs + (1+log(n)) * npars
    names(devs) <- models
    indx <- switch(match.arg(orderby),
                AIC = sort.list(criteria[1,]),
                BIC = sort.list(criteria[2,]),
                SC = sort.list(criteria[3,]))
    estimate <- estimate[, indx]
    std.errors <- std.errors[, indx]
    criteria <- criteria[, indx]
    devs <- devs[indx]
    dm <- dm[, indx]
    list(estimate = estimate, std.err = std.errors, deviance = devs,
         criteria = criteria, dep.summary = dm)
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

na.vals <- function(x) {
    x <- cbind(is.na(x[,1]),2*is.na(x[,2]))
    drop(x %*% c(1,1))
}


