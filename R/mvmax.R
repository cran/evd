"rmvlog"<-
# Uses Algorithm 2.1 in Stephenson(2002)
function(n, dep, d = 2, mar = c(0,1,0))
{
    if(length(dep) != 1 || mode(dep) != "numeric" || dep <= 0 ||
       dep > 1) stop("invalid argument for `dep'")
    #ptm <- proc.time()
    sim <- .C("rmvlog_tawn",
              as.integer(n), as.integer(d), as.double(dep),
              sim = double(d*n), PACKAGE = "evd")$sim
    #return(proc.time()-ptm)
    mtransform(matrix(1/sim, ncol=d, byrow=TRUE), mar, inv = TRUE)
}

"rmvalog"<-
# Uses Algorithm 2.2 in Stephenson(2002)
function(n, dep, asy, d = 2, mar = c(0,1,0))
{
    nb <- 2^d-1
    dep <- rep(dep, length.out = nb-d)
    if(any(dep <= 0) || any(dep > 1))
        stop("invalid argument for `dep'")
    # Poor Interface
    subsets <- function(d) {
        x <- 1:d
        k <- NULL
        for(m in x) 
            k <- rbind(cbind(TRUE, k), cbind(FALSE, k))
        pset <- apply(k, 1, function(z) x[z])
        pset[sort.list(unlist(lapply(pset,length)))[-1]] 
    }
    tasy <- function(theta) {
        b <- subsets(d)
        trans <- matrix(0,nrow=nb,ncol=d)
        for(i in 1:nb) {
            if(theta == "equal")
                trans[i,(1:d %in% b[[i]])] <- 1/2^(d-1)
            else
            trans[i,(1:d %in% b[[i]])] <- theta[[i]]
        }
        trans
    }
    if(mode(asy) == "list" && length(asy) == nb || asy == "equal")
        asy <- tasy(asy)
    if(!is.matrix(asy) || mode(asy) != "numeric" || any(dim(asy) != c(nb,d)))
        stop("`asy' is not of the correct form")
    if(min(asy) < 0 || max(asy) > 1)
       stop("`asy' must contain parameters in [0,1]")
    if(any(apply(asy,2,sum) != 1) || any(asy[c(rep(FALSE,d),dep==1),]!=0) ||
       any(apply(asy[-(1:d),,drop=FALSE],1,function(x) sum(x!=0)) == 1))
       stop("`asy' does not satisfy the appropriate constraints")
    dep <- c(rep(1,d),dep)
    #ptm <- proc.time()
    sim <- .C("rmvalog_tawn",
              as.integer(n), as.integer(d), as.integer(nb), as.double(dep),
              as.double(t(asy)), sim = double(n*d), PACKAGE = "evd")$sim
    #return(proc.time()-ptm)
    mtransform(matrix(1/sim, ncol=d, byrow=TRUE), mar, inv = TRUE)
}

"pmvlog"<- 
function(q, dep, d = 2, mar = c(0,1,0))
{
    if(length(dep) != 1 || mode(dep) != "numeric" || dep <= 0 ||
        dep > 1) stop("invalid argument for `dep'")
    if(is.null(dim(q))) dim(q) <- c(1,d)
    if(ncol(q) != d) stop("`q' and `d' are not compatible")
    q <- mtransform(q, mar)
    exp(-apply(q^(1/dep),1,sum)^dep)
}

"pmvalog"<-
function(q, dep, asy, d = 2, mar = c(0,1,0))
{
    nb <- 2^d-1
    dep <- rep(dep, length.out = nb-d)
    if(any(dep <= 0) || any(dep > 1))
        stop("invalid argument for `dep'")
    # Poor Interface
    subsets <- function(d) {
        x <- 1:d
        k <- NULL
        for(m in x) 
            k <- rbind(cbind(TRUE, k), cbind(FALSE, k))
        pset <- apply(k, 1, function(z) x[z])
        pset[sort.list(unlist(lapply(pset,length)))[-1]] 
    }
    tasy <- function(theta) {
        b <- subsets(d)
        trans <- matrix(0,nrow=nb,ncol=d)
        for(i in 1:nb) {
            if(theta == "equal")
                trans[i,(1:d %in% b[[i]])] <- 1/2^(d-1)
            else
            trans[i,(1:d %in% b[[i]])] <- theta[[i]]
        }
        trans
    }
    if(mode(asy) == "list" && length(asy) == nb || asy == "equal")
        asy <- tasy(asy)
    if(!is.matrix(asy) || mode(asy) != "numeric" || any(dim(asy) != c(nb,d)))
        stop("`asy' is not of the correct form")
    if(min(asy) < 0 || max(asy) > 1)
       stop("`asy' must contain parameters in [0,1]")
    if(any(apply(asy,2,sum) != 1) || any(asy[c(rep(FALSE,d),dep==1),]!=0) ||
       any(apply(asy[-(1:d),,drop=FALSE],1,function(x) sum(x!=0)) == 1))
       stop("`asy' does not satisfy the appropriate constraints")
    dep <- c(rep(1,d),dep)
    if(is.null(dim(q))) dim(q) <- c(1,d)
    if(ncol(q) != d) stop("`q' and `d' are not compatible")
    q <- mtransform(q, mar)
    inner <- function(par)
        apply((rep(par[1:d], rep(nrow(q),d))*q)^(1/par[d+1]), 1, sum)^par[d+1]
    comps <- apply(cbind(asy,dep),1,inner)
    if(is.null(dim(comps))) dim(comps) <- c(1,nb)
    exp(-apply(comps,1,sum))
}













