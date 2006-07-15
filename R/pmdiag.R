"mrlplot"<-
function(data, tlim, pscale = FALSE, nt = max(100, length(data)), lty =
    c(2,1,2), col = 1, conf = 0.95, main = "Mean Residual Life Plot",
    xlab = "Threshold", ylab = "Mean Excess", ...)
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
    if(pscale) { 
      tlim[1] <- mean(data <= tlim[1], na.rm = TRUE)
      tlim[2] <- mean(data <= tlim[2], na.rm = TRUE)
      pvec <- seq(tlim[1], tlim[2], length = nt)
      u <- quantile(data, probs = pvec, na.rm = TRUE)
    }
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
    if(pscale) {
      u <- pvec
      if(missing(xlab)) xlab <- "Threshold probability"
    }
    matplot(u, x, type = "l", lty = lty, col = col, main = main,
            xlab = xlab, ylab = ylab, ...)
    invisible(list(x = u, y = x))
}

"tcplot"<-
function(data, tlim, model = c("gpd", "pp"), pscale = FALSE, cmax = FALSE,
    r = 1, ulow = -Inf, rlow = 1, nt = 25, which = 1:npar, conf = 0.95,
    lty = 1, lwd = 1, type = "b", cilty = 1, vci = TRUE, xlab, xlim, ylabs,
    ylims, ask = nb.fig < length(which) && dev.interactive(), ...)
{
    model <- match.arg(model)
    u <- seq(tlim[1], tlim[2], length = nt)
    if(pscale) { 
      tlim[1] <- mean(data <= tlim[1], na.rm = TRUE)
      tlim[2] <- mean(data <= tlim[2], na.rm = TRUE)
      pvec <- seq(tlim[1], tlim[2], length = nt)
      u <- quantile(data, probs = pvec, na.rm = TRUE)
    }
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
   if(pscale) u <- pvec
   if(missing(xlim)) xlim <- tlim
   if(missing(xlab)) {
     xlab <- "Threshold"
     if(pscale) xlab <- "Threshold probability"
   }
   if(model == "pp") {
     ylab <- c("Location","Scale","Shape")
     if(!missing(ylabs)) ylab[show] <- ylabs
     ylim <- rbind(range(locs), range(scls), range(shps))
     if(!missing(ylims)) ylim[show,] <- ylims
     if(show[1]) {
       matplot(u, locs, type = "n", xlab = xlab, ylab = ylab[1],
            xlim = xlim, ylim = ylim[1,])
       lines(u, locs[,2], lty = lty, lwd = lwd, type = type)
       if(vci) segments(u, locs[,1], u, locs[,3], lty = cilty)
       else {
         lines(u, locs[,1], lty = cilty)
         lines(u, locs[,3], lty = cilty)
       }
     }
     if(show[2]) {
       matplot(u, scls, type = "n", xlab = xlab, ylab = ylab[2],
            xlim = xlim, ylim = ylim[2,])
       lines(u, scls[,2], lty = lty, lwd = lwd, type = type)
       if(vci) segments(u, scls[,1], u, scls[,3], lty = cilty)
       else {
         lines(u, scls[,1], lty = cilty)
         lines(u, scls[,3], lty = cilty)
       }
     }
     if(show[3]) {
       matplot(u, shps, type = "n", xlab = xlab, ylab = ylab[3],
            xlim = xlim, ylim = ylim[3,])
       lines(u, shps[,2], lty = lty, lwd = lwd, type = type)
       if(vci) segments(u, shps[,1], u, shps[,3], lty = cilty)
       else {
         lines(u, shps[,1], lty = cilty)
         lines(u, shps[,3], lty = cilty)
       }
     }
     rtlist <- list(locs = locs, scales = scls, shapes = shps) 
   }
   if(model == "gpd") {
     ylab <- c("Modified Scale","Shape")
     if(!missing(ylabs)) ylab[show] <- ylabs
     ylim <- rbind(range(scls), range(shps))
     if(!missing(ylims)) ylim[show,] <- ylims
     if(show[1]) {
       matplot(u, scls, type = "n", xlab = xlab, ylab = ylab[1],
            xlim = xlim, ylim = ylim[1,])
       lines(u, scls[,2], lty = lty, lwd = lwd, type = type)
       if(vci) segments(u, scls[,1], u, scls[,3], lty = cilty)
       else {
         lines(u, scls[,1], lty = cilty)
         lines(u, scls[,3], lty = cilty)
       }
     }
     if(show[2]) {
       matplot(u, shps, type = "n", xlab = xlab, ylab = ylab[2],
            xlim = xlim, ylim = ylim[2,])
       lines(u, shps[,2], lty = lty, lwd = lwd, type = type)
       if(vci) segments(u, shps[,1], u, shps[,3], lty = cilty)
       else {
         lines(u, shps[,1], lty = cilty)
         lines(u, shps[,3], lty = cilty)
       }
     }
     rtlist <- list(scales = scls, shapes = shps) 
   } 
   invisible(rtlist)
}

chiplot <- function(data, nq = 100, qlim = NULL, which = 1:2, conf = 0.95, boot = FALSE, spcases = FALSE, lty = 1, cilty = 2, col = 1, cicol = 1, xlim = c(0,1), ylim1 = NULL, ylim2 = c(-1,1), main1 = "Chi Plot", main2 = "Chi Bar Plot", xlab = "Quantile", ylab1 = "Chi", ylab2 = "Chi Bar", ask = nb.fig < length(which) && dev.interactive(), ...)
{
    data <- na.omit(data)
    n <- nrow(data)
    data <- cbind(rank(data[, 1])/(n + 1), rank(data[, 2])/(n + 1))
    rowmax <- apply(data, 1, max)
    rowmin <- apply(data, 1, min)
    eps <- .Machine$double.eps^0.5
    qlim2 <- c(min(rowmax) + eps, max(rowmin) - eps)
    if(!is.null(qlim)) {
      if(qlim[1] < qlim2[1]) stop("lower quantile limit is too low")
      if(qlim[2] > qlim2[2]) stop("upper quantile limit is too high")
      if(qlim[1] > qlim[2]) stop("lower quantile limit is less than upper quantile limit")
    } else qlim <- qlim2
    u <- seq(qlim[1], qlim[2], length = nq)

    cu <- cbaru <- numeric(nq)
    for(i in 1:nq) cu[i] <- mean(rowmax < u[i])
    for(i in 1:nq) cbaru[i] <- mean(rowmin > u[i])
    chiu <- 2 - log(cu)/log(u)
    chibaru <- (2 * log(1 - u))/log(cbaru) - 1

    if(!boot) {
      cnst <- qnorm((1 + conf)/2)
      varchi <- ((1/log(u)^2 * 1)/cu^2 * cu * (1 - cu))/n
      varchi <- cnst * sqrt(varchi)
      varchibar <- (((4 * log(1 - u)^2)/(log(cbaru)^4 * cbaru^2)) * cbaru * (
    	1 - cbaru))/n
      varchibar <- cnst * sqrt(varchibar)
      chiu <- cbind(chilow = chiu-varchi, chi = chiu, chiupp = chiu+varchi) 
      chibaru <- cbind(chiblow = chibaru-varchibar, chib = chibaru, chibupp =
        chibaru+varchibar)
    }
    else {
      cui <- cbarui <- matrix(0, nrow = nq, ncol = 100)
      for(i in 1:100) {
        datai <- data[sample(n, replace = TRUE),]
        rowmaxi <- apply(datai, 1, max)
        rowmini <- apply(datai, 1, min)
        for(j in 1:nq) cui[j,i] <- mean(rowmaxi < u[j])
        for(j in 1:nq) cbarui[j,i] <- mean(rowmini > u[j])
      }
      chiucon <- 2 - log(cui)/log(u)
      chibarucon <- (2 * log(1 - u))/log(cbarui) - 1
      chiucon <- apply(chiucon, 1, quantile, probs = c(0.025, 0.975))
      chibarucon <- apply(chibarucon, 1, quantile, probs = c(0.025, 0.975))
      chiu <- cbind(chilow = chiucon[1,], chi = chiu, chiupp = chiucon[2,]) 
      chibaru <- cbind(chiblow = chibarucon[1,], chib = chibaru, chibupp =
        chibarucon[2,])
    }

    show <- logical(2)
    show[which] <- TRUE
    lty <- c(cilty, lty, cilty)
    col <- c(cicol, col, cicol)
    nb.fig <- prod(par("mfcol"))
    if (ask) {
      op <- par(ask = TRUE)
      on.exit(par(op))
    }
    
    if(is.null(ylim1)) ylim1 <- c(min(c(chiu, 0)), 1)
    if(show[1]) {
      matplot(u, chiu, type = "l", lty = lty, col = col, xlim = xlim, ylim = ylim1,
              main = main1, xlab = xlab, ylab = ylab1, ...)
      if(spcases) {
        segments(0,0,1,0, lty = 5, col = "grey")
        segments(0,1,1,1, lty = 5, col = "grey")
        lines(u, 2-log(pmax(2*u-1,0))/log(u), lty = 5, col = "grey")
      }
    }

    if(show[2]) {
      matplot(u, chibaru, type = "l", lty = lty, col = col, xlim = xlim, ylim = ylim2,
              main = main2, xlab = xlab, ylab = ylab2, ...)
      if(spcases) {
        segments(0,0,1,0, lty = 5, col = "grey")
        segments(0,1,1,1, lty = 5, col = "grey")
        lines(u, 2*log(1-u)/log(1-2*u+pmax(2*u-1,0)) - 1, lty = 5, col = "grey")
      }
    }

    plvals <- list(quantile = u, chi = chiu, chibar = chibaru)
    if(!show[1]) plvals$chi <- NULL
    if(!show[2]) plvals$chib <- NULL
    invisible(plvals)    
}



