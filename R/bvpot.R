
################## bivariate threshold fitting routines #################

fbvpot <- function(x, threshold, model = c("log", "bilog", "alog", "neglog", "negbilog", "aneglog", "ct"), likelihood = c("censored","poisson"), start, ..., sym = FALSE, std.err = TRUE, dsm = TRUE, corr = FALSE, method = "BFGS", warn.inf = TRUE) {
  call <- match.call()
  likelihood <- match.arg(likelihood)
  ft <- switch(likelihood,
    censored = fbvcpot(x = x, u = threshold, model = model, start = start, ...,
      sym = sym, std.err = std.err, dsm = dsm, corr = corr,
      method = method, warn.inf = warn.inf),
    poisson = fbvppot(x = x, u = threshold, model = model, start = start, ...,
      sym = sym, std.err = std.err, dsm = dsm, corr = corr,
      method = method, warn.inf = warn.inf))
  structure(c(ft, call = call), class = c("bvpot", "evd"))
}

fbvcpot <- function(x, u, model = c("log", "bilog", "alog", "neglog", "negbilog", "aneglog", "ct"), start, ..., sym = FALSE, std.err = TRUE, dsm = TRUE, corr = FALSE, method = "BFGS", warn.inf = TRUE) {
  model <- match.arg(model)
  if(sym && !(model %in% c("alog","aneglog","ct")))
    warning("Argument `sym' was ignored")
  switch(model,
    log = fbvclog(x = x, u = u, start = start, ..., std.err = std.err,
      dsm = dsm, corr = corr, method = method, warn.inf = warn.inf),
    bilog = fbvcbilog(x = x, u = u, start = start, ..., std.err = std.err,
      dsm = dsm, corr = corr, method = method, warn.inf = warn.inf),
    alog = fbvcalog(x = x, u = u, start = start, ..., sym = sym,
      std.err = std.err, dsm = dsm, corr = corr, method = method, warn.inf =
      warn.inf),
    neglog = fbvcneglog(x = x, u = u, start = start, ..., std.err = std.err,
      dsm = dsm, corr = corr, method = method, warn.inf = warn.inf),
    negbilog = fbvcnegbilog(x = x, u = u, start = start, ..., std.err =
      std.err, dsm = dsm, corr = corr, method = method, warn.inf = warn.inf),
    aneglog = fbvcaneglog(x = x, u = u, start = start, ..., sym =
      sym, std.err = std.err, dsm = dsm, corr = corr, method = method,
      warn.inf = warn.inf),
    ct = fbvcct(x = x, u = u, start = start, ..., sym = sym,
      std.err = std.err, dsm = dsm, corr = corr, method = method, warn.inf =
      warn.inf))
}

fbvppot <- function(x, u, model = c("log", "bilog", "alog", "neglog", "negbilog", "aneglog", "ct"), start, ..., sym = FALSE, std.err = TRUE, dsm = TRUE, corr = FALSE, method = "BFGS", warn.inf = TRUE) {
  model <- match.arg(model)
  if(model %in% c("alog","aneglog"))
    stop("This model is not implemented for poisson likelihood")
  if(sym && (model != "ct"))
    warning("Argument `sym' was ignored")
  switch(model,
    log = fbvplog(x = x, u = u, start = start, ..., std.err = std.err,
      dsm = dsm, corr = corr, method = method, warn.inf = warn.inf),
    bilog = fbvpbilog(x = x, u = u, start = start, ..., std.err = std.err,
      dsm = dsm, corr = corr, method = method, warn.inf = warn.inf),
    neglog = fbvpneglog(x = x, u = u, start = start, ..., std.err = std.err,
      dsm = dsm, corr = corr, method = method, warn.inf = warn.inf),
    negbilog = fbvpnegbilog(x = x, u = u, start = start, ..., std.err =
      std.err, dsm = dsm, corr = corr, method = method, warn.inf = warn.inf),
    ct = fbvpct(x = x, u = u, start = start, ..., sym = sym,
      std.err = std.err, dsm = dsm, corr = corr, method = method, warn.inf =
      warn.inf))
}

################## censored likelihood fitting routines #################

fbvclog <- function(x, u, start, ..., std.err = TRUE, dsm = TRUE, corr = FALSE, method = "BFGS", warn.inf = TRUE) {

  nllbvclog <- function(scale1, shape1, scale2, shape2, dep) { 
    .C("nllbvclog", spx$x1, spx$x2, spx$nn, spx$n, spx$thdi,
      spx$lambda, dep, scale1, shape1, scale2, shape2, dns = double(1),
      PACKAGE = "evd")$dns
  }
  param <- c("scale1", "shape1", "scale2", "shape2", "dep")
  nmdots <- names(list(...))
  start <- bvstart.vals(x, start, NULL, NULL, nmdots, param, NULL,
    NULL, model = "log", sym = sym, obj = "bvpot", u = u)
  spx <- sep.bvdata(x, obj = "bvpot", u = u)  
  nm <- names(start)
  fixed.param <- list(...)[nmdots %in% param]
  if(any(!(param %in% c(nm,names(fixed.param)))))
    stop("unspecified parameters")
  l <- length(nm)
  f <- formals(nllbvclog)
  names(f) <- param
  m <- match(nm, param)
  if(any(is.na(m))) 
    stop("`start' specifies unknown arguments")     
  formals(nllbvclog) <- c(f[m], f[-m])
  nll <- function(p, ...) nllbvclog(p, ...)
  if(l > 1) {
    body(nll) <- parse(text = paste("nllbvclog(", paste("p[", 1:l, "]", collapse=", "), ", ...)"))
  }
  start.arg <- c(list(p = unlist(start)), fixed.param)
  if(warn.inf && do.call("nll", start.arg) == 1e+06)
    warning("negative log-likelihood is infinite at starting values")
  opt <- optim(start, nll, hessian = std.err, ..., method = method)
  bvtpost.optim(x, u, opt, nm, fixed.param, std.err, dsm, corr, spx$nat, FALSE, model = "log")
}

fbvcbilog <- function(x, u, start, ..., std.err = TRUE, dsm = TRUE, corr = FALSE, method = "BFGS", warn.inf = TRUE) {

  nllbvcbilog <- function(scale1, shape1, scale2, shape2, alpha, beta) {
    .C("nllbvcbilog", spx$x1, spx$x2, spx$nn, spx$n, spx$thdi, spx$lambda,
      alpha, beta, scale1, shape1, scale2, shape2, dns = double(1),
      PACKAGE = "evd")$dns
  }
  param <- c("scale1", "shape1", "scale2", "shape2", "alpha", "beta")
  nmdots <- names(list(...))
  start <- bvstart.vals(x, start, NULL, NULL, nmdots, param, NULL,
    NULL, model = "bilog", sym = sym, obj = "bvpot", u = u)
  spx <- sep.bvdata(x, obj = "bvpot", u = u) 
  nm <- names(start)
  fixed.param <- list(...)[nmdots %in% param]
  if(any(!(param %in% c(nm,names(fixed.param)))))
    stop("unspecified parameters")
  l <- length(nm)
  f <- formals(nllbvcbilog)
  names(f) <- param
  m <- match(nm, param)
  if(any(is.na(m))) 
    stop("`start' specifies unknown arguments")   
  formals(nllbvcbilog) <- c(f[m], f[-m])
  nll <- function(p, ...) nllbvcbilog(p, ...)
  if(l > 1) {
    body(nll) <- parse(text = paste("nllbvcbilog(", paste("p[", 1:l,
      "]", collapse=", "), ", ...)"))
  }
  start.arg <- c(list(p = unlist(start)), fixed.param)
  if(warn.inf && do.call("nll", start.arg) == 1e+06)
    warning("negative log-likelihood is infinite at starting values")
  opt <- optim(start, nll, hessian = std.err, ..., method = method)
  bvtpost.optim(x, u, opt, nm, fixed.param, std.err, dsm, corr, spx$nat, FALSE, model = "bilog")
}

fbvcalog <- function(x, u, start, ..., sym = FALSE, std.err = TRUE, dsm = TRUE, corr = FALSE, method = "BFGS", warn.inf = TRUE) {

  nllbvcalog.sym <- function(scale1, shape1, scale2, shape2, asy1, dep) {
    nllbvcalog(scale1, shape1, scale2, shape2, asy1, asy1, dep)
  }

  nllbvcalog <- function(scale1, shape1, scale2, shape2, asy1, asy2, dep) {
    .C("nllbvcalog", spx$x1, spx$x2, spx$nn, spx$n, spx$thdi, spx$lambda,
       dep, asy1, asy2, scale1, shape1, scale2, shape2, dns = double(1),
       PACKAGE = "evd")$dns
  }
  param <- c("scale1", "shape1", "scale2", "shape2", "asy1", "dep")
  if(!sym) param[6:7] <- c("asy2", "dep")
  nmdots <- names(list(...))
  start <- bvstart.vals(x, start, NULL, NULL, nmdots, param, NULL,
    NULL, model = "alog", sym = sym, obj = "bvpot", u = u)
  spx <- sep.bvdata(x, obj = "bvpot", u = u)
  nm <- names(start)
  fixed.param <- list(...)[nmdots %in% param]
  if(any(!(param %in% c(nm,names(fixed.param)))))
    stop("unspecified parameters")
  l <- length(nm)
  if(sym) f <- formals(nllbvcalog.sym)
  else f <- formals(nllbvcalog)
  names(f) <- param
  m <- match(nm, param)
  if(any(is.na(m))) 
    stop("`start' specifies unknown arguments") 
  if(sym) {
    formals(nllbvcalog.sym) <- c(f[m], f[-m])
    nll <- function(p, ...) nllbvcalog.sym(p, ...)
    if(l > 1) {
      body(nll) <- parse(text = paste("nllbvcalog.sym(", paste("p[", 1:l, "]", collapse=", "), ", ...)"))
    }
  } else {
    formals(nllbvcalog) <- c(f[m], f[-m])
    nll <- function(p, ...) nllbvcalog(p, ...)
    if(l > 1) {
      body(nll) <- parse(text = paste("nllbvcalog(", paste("p[", 1:l, "]", collapse=", "), ", ...)"))
    }
  }
  start.arg <- c(list(p = unlist(start)), fixed.param)
  if(warn.inf && do.call("nll", start.arg) == 1e+06)
    warning("negative log-likelihood is infinite at starting values")
  opt <- optim(start, nll, hessian = std.err, ..., method = method)
  if(sym) {
    asy2 <- opt$par["asy1"]
    names(asy2) <- NULL
    if(is.na(asy2)) asy2 <- fixed.param[["asy1"]]
    fixed.param <- c(fixed.param, list(asy2 = asy2))
  }
  bvtpost.optim(x, u, opt, nm, fixed.param, std.err, dsm, corr,
    spx$nat, sym, model = "alog")
}

fbvcneglog <- function(x, u, start, ..., std.err = TRUE, dsm = TRUE, corr = FALSE, method = "BFGS", warn.inf = TRUE) {

  nllbvcneglog <- function(scale1, shape1, scale2, shape2, dep) {
    .C("nllbvcneglog", spx$x1, spx$x2, spx$nn, spx$n, spx$thdi, spx$lambda,
      dep, scale1, shape1, scale2, shape2, dns = double(1),
      PACKAGE = "evd")$dns
  }
  param <- c("scale1", "shape1", "scale2", "shape2", "dep")
  nmdots <- names(list(...))
  start <- bvstart.vals(x, start, NULL, NULL, nmdots, param, NULL,
    NULL, model = "neglog", sym = sym, obj = "bvpot", u = u)
  spx <- sep.bvdata(x, obj = "bvpot", u = u) 
  nm <- names(start)
  fixed.param <- list(...)[nmdots %in% param]
  if(any(!(param %in% c(nm,names(fixed.param)))))
    stop("unspecified parameters")
  l <- length(nm)
  f <- formals(nllbvcneglog)
  names(f) <- param
  m <- match(nm, param)
  if(any(is.na(m))) 
    stop("`start' specifies unknown arguments") 
  formals(nllbvcneglog) <- c(f[m], f[-m])
  nll <- function(p, ...) nllbvcneglog(p, ...)
  if(l > 1) {
    body(nll) <- parse(text = paste("nllbvcneglog(", paste("p[", 1:l, "]", collapse=", "), ", ...)"))
  }
  start.arg <- c(list(p = unlist(start)), fixed.param)
  if(warn.inf && do.call("nll", start.arg) == 1e+06)
    warning("negative log-likelihood is infinite at starting values")
  opt <- optim(start, nll, hessian = std.err, ..., method = method)
  bvtpost.optim(x, u, opt, nm, fixed.param, std.err, dsm, corr, spx$nat, FALSE, model = "neglog")
}

fbvcnegbilog <- function(x, u, start, ..., std.err = TRUE, dsm = TRUE, corr = FALSE, method = "BFGS", warn.inf = TRUE) {

  nllbvcnegbilog <- function(scale1, shape1, scale2, shape2, alpha, beta) {
    .C("nllbvcnegbilog", spx$x1, spx$x2, spx$nn, spx$n, spx$thdi, spx$lambda,
       alpha, beta, scale1, shape1, scale2, shape2, dns = double(1),
       PACKAGE = "evd")$dns
  }
  param <- c("scale1", "shape1", "scale2", "shape2", "alpha", "beta")
  nmdots <- names(list(...))
  start <- bvstart.vals(x, start, NULL, NULL, nmdots, param, NULL,
    NULL, model = "negbilog", sym = sym, obj = "bvpot", u = u)
  spx <- sep.bvdata(x, obj = "bvpot", u = u) 
  nm <- names(start)
  fixed.param <- list(...)[nmdots %in% param]
  if(any(!(param %in% c(nm,names(fixed.param)))))
    stop("unspecified parameters")
  l <- length(nm)
  f <- formals(nllbvcnegbilog)
  names(f) <- param
  m <- match(nm, param)
  if(any(is.na(m))) 
    stop("`start' specifies unknown arguments") 
  formals(nllbvcnegbilog) <- c(f[m], f[-m])
  nll <- function(p, ...) nllbvcnegbilog(p, ...)
  if(l > 1) {
    body(nll) <- parse(text = paste("nllbvcnegbilog(", paste("p[", 1:l,
      "]", collapse=", "), ", ...)"))
  }
  start.arg <- c(list(p = unlist(start)), fixed.param)
  if(warn.inf && do.call("nll", start.arg) == 1e+06)
    warning("negative log-likelihood is infinite at starting values")
  opt <- optim(start, nll, hessian = std.err, ..., method = method)
  bvtpost.optim(x, u, opt, nm, fixed.param, std.err, dsm, corr, spx$nat, FALSE, model = "negbilog")
}

fbvcaneglog <- function(x, u, start, ..., sym = FALSE, std.err = TRUE, dsm = TRUE, corr = FALSE, method = "BFGS", warn.inf = TRUE) {

  nllbvcaneglog.sym <- function(scale1, shape1, scale2, shape2, asy1, dep) {
    nllbvcaneglog(scale1, shape1, scale2, shape2, asy1, asy1, dep)
  }

  nllbvcaneglog <- function(scale1, shape1, scale2, shape2, asy1, asy2, dep) {
    .C("nllbvcaneglog", spx$x1, spx$x2, spx$nn, spx$n, spx$thdi, spx$lambda,
       dep, asy1, asy2, scale1, shape1, scale2, shape2, dns = double(1),
       PACKAGE = "evd")$dns
  }
  param <- c("scale1", "shape1", "scale2", "shape2", "asy1", "dep")
  if(!sym) param[6:7] <- c("asy2", "dep")
  nmdots <- names(list(...))
  start <- bvstart.vals(x, start, NULL, NULL, nmdots, param, NULL,
    NULL, model = "aneglog", sym = sym, obj = "bvpot", u = u)
  spx <- sep.bvdata(x, obj = "bvpot", u = u) 
  nm <- names(start)
  fixed.param <- list(...)[nmdots %in% param]
  if(any(!(param %in% c(nm,names(fixed.param)))))
    stop("unspecified parameters")
  l <- length(nm)
  if(sym) f <- formals(nllbvcaneglog.sym)
  else f <- formals(nllbvcaneglog)
  names(f) <- param
  m <- match(nm, param)
  if(any(is.na(m))) 
    stop("`start' specifies unknown arguments") 
  if(sym) {
    formals(nllbvcaneglog.sym) <- c(f[m], f[-m])
    nll <- function(p, ...) nllbvcaneglog.sym(p, ...)
    if(l > 1) {
      body(nll) <- parse(text = paste("nllbvcaneglog.sym(", paste("p[", 1:l, "]", collapse=", "), ", ...)"))
    }
  } else {
    formals(nllbvcaneglog) <- c(f[m], f[-m])
    nll <- function(p, ...) nllbvcaneglog(p, ...)
    if(l > 1) {
      body(nll) <- parse(text = paste("nllbvcaneglog(", paste("p[", 1:l, "]", collapse=", "), ", ...)"))
    }
  }
  start.arg <- c(list(p = unlist(start)), fixed.param)
  if(warn.inf && do.call("nll", start.arg) == 1e+06)
    warning("negative log-likelihood is infinite at starting values")
  opt <- optim(start, nll, hessian = std.err, ..., method = method)
  if(sym) {
    asy2 <- opt$par["asy1"]
    names(asy2) <- NULL
    if(is.na(asy2)) asy2 <- fixed.param[["asy1"]]
    fixed.param <- c(fixed.param, list(asy2 = asy2))
  }
  bvtpost.optim(x, u, opt, nm, fixed.param, std.err, dsm, corr, spx$nat, sym, model = "aneglog")
}

fbvcct <- function(x, u, start, ..., sym = FALSE, std.err = TRUE, dsm = TRUE, corr = FALSE, method = "BFGS", warn.inf = TRUE) {

  nllbvcct.sym <- function(scale1, shape1, scale2, shape2, alpha) {
    nllbvcct(scale1, shape1, scale2, shape2, alpha, alpha)
  }

  nllbvcct <- function(scale1, shape1, scale2, shape2, alpha, beta) {
    .C("nllbvcct", spx$x1, spx$x2, spx$nn, spx$n, spx$thdi, spx$lambda,
      alpha, beta, scale1, shape1, scale2, shape2, dns = double(1),
      PACKAGE = "evd")$dns
  }
  param <- c("scale1", "shape1", "scale2", "shape2", "alpha")
  if(!sym) param <- c(param, "beta")
  nmdots <- names(list(...))
  start <- bvstart.vals(x, start, NULL, NULL, nmdots, param, NULL,
    NULL, model = "ct", sym = sym, obj = "bvpot", u = u)
  spx <- sep.bvdata(x, obj = "bvpot", u = u) 
  nm <- names(start)
  fixed.param <- list(...)[nmdots %in% param]
  if(any(!(param %in% c(nm,names(fixed.param)))))
    stop("unspecified parameters")
  l <- length(nm)
  if(sym) f <- formals(nllbvcct.sym)
  else f <- formals(nllbvcct)
  names(f) <- param
  m <- match(nm, param)
  if(any(is.na(m))) 
    stop("`start' specifies unknown arguments") 
  if(sym) {
    formals(nllbvcct.sym) <- c(f[m], f[-m])
    nll <- function(p, ...) nllbvcct.sym(p, ...)
    if(l > 1) {
      body(nll) <- parse(text = paste("nllbvcct.sym(", paste("p[", 1:l, "]", collapse=", "), ", ...)"))
    }
  } else {
    formals(nllbvcct) <- c(f[m], f[-m])
    nll <- function(p, ...) nllbvcct(p, ...)
    if(l > 1) {
      body(nll) <- parse(text = paste("nllbvcct(", paste("p[", 1:l, "]", collapse=", "), ", ...)"))
    }
  }
  start.arg <- c(list(p = unlist(start)), fixed.param)
  if(warn.inf && do.call("nll", start.arg) == 1e+06)
    warning("negative log-likelihood is infinite at starting values")
  opt <- optim(start, nll, hessian = std.err, ..., method = method)
  if(sym) {
    beta <- opt$par["alpha"]
    names(beta) <- NULL
    if(is.na(beta)) beta <- fixed.param[["alpha"]]
    fixed.param <- c(fixed.param, list(beta = beta))
  }
  bvtpost.optim(x, u, opt, nm, fixed.param, std.err, dsm, corr, spx$nat, sym, model = "ct")
}


################## Poisson likelihood fitting routines ##################

fbvplog <- function(x, u, start, ..., std.err = TRUE, dsm = TRUE, corr = FALSE, method = "BFGS", warn.inf = TRUE) {

  nllbvplog <- function(scale1, shape1, scale2, shape2, dep) {
    .C("nllbvplog", spx$x1, spx$x2, spx$nn, spx$n, spx$thdi, spx$r1, spx$r2,
      spx$lambda, dep, scale1, shape1, scale2, shape2, dns = double(1),
      PACKAGE = "evd")$dns
  }
  param <- c("scale1", "shape1", "scale2", "shape2", "dep")
  nmdots <- names(list(...))
  start <- bvstart.vals(x, start, NULL, NULL, nmdots, param, NULL,
    NULL, model = "log", sym = sym, obj = "bvpot", u = u)
  spx <- sep.bvdata(x, obj = "bvpot", u = u, censored = FALSE) 
  nm <- names(start)
  fixed.param <- list(...)[nmdots %in% param]
  if(any(!(param %in% c(nm,names(fixed.param)))))
    stop("unspecified parameters")
  l <- length(nm)
  f <- formals(nllbvplog)
  names(f) <- param
  m <- match(nm, param)
  if(any(is.na(m))) 
    stop("`start' specifies unknown arguments")   
  formals(nllbvplog) <- c(f[m], f[-m])
  nll <- function(p, ...) nllbvplog(p, ...)
  if(l > 1) {
    body(nll) <- parse(text = paste("nllbvplog(", paste("p[", 1:l, "]", collapse=", "), ", ...)"))
  }
  start.arg <- c(list(p = unlist(start)), fixed.param)
  if(warn.inf && do.call("nll", start.arg) == 1e+06)
    warning("negative log-likelihood is infinite at starting values")
  opt <- optim(start, nll, hessian = std.err, ..., method = method)
  bvtpost.optim(x, u, opt, nm, fixed.param, std.err, dsm, corr,
    spx$nat, FALSE, model = "log")
}

fbvpneglog <- function(x, u, start, ..., std.err = TRUE, dsm = TRUE, corr = FALSE, method = "BFGS", warn.inf = TRUE) {

  nllbvpneglog <- function(scale1, shape1, scale2, shape2, dep) {
    .C("nllbvpneglog", spx$x1, spx$x2, spx$nn, spx$n, spx$thdi, spx$r1,
      spx$r2, spx$lambda, dep, scale1, shape1, scale2, shape2,
      dns = double(1), PACKAGE = "evd")$dns    
  }
  param <- c("scale1", "shape1", "scale2", "shape2", "dep")
  nmdots <- names(list(...))
  start <- bvstart.vals(x, start, NULL, NULL, nmdots, param, NULL,
    NULL, model = "neglog", sym = sym, obj = "bvpot", u = u)
  spx <- sep.bvdata(x, obj = "bvpot", u = u, censored = FALSE)
  nm <- names(start)
  fixed.param <- list(...)[nmdots %in% param]
  if(any(!(param %in% c(nm,names(fixed.param)))))
    stop("unspecified parameters")
  l <- length(nm)
  f <- formals(nllbvpneglog)
  names(f) <- param
  m <- match(nm, param)
  if(any(is.na(m))) 
    stop("`start' specifies unknown arguments") 
  formals(nllbvpneglog) <- c(f[m], f[-m])
  nll <- function(p, ...) nllbvpneglog(p, ...)
  if(l > 1) {
    body(nll) <- parse(text = paste("nllbvpneglog(", paste("p[", 1:l, "]", collapse=", "), ", ...)"))
  }
  start.arg <- c(list(p = unlist(start)), fixed.param)
  if(warn.inf && do.call("nll", start.arg) == 1e+06)
    warning("negative log-likelihood is infinite at starting values")
  opt <- optim(start, nll, hessian = std.err, ..., method = method)
  bvtpost.optim(x, u, opt, nm, fixed.param, std.err, dsm, corr, spx$nat, FALSE, model = "neglog")
}

fbvpct <- function(x, u, start, ..., sym = FALSE, std.err = TRUE, dsm = TRUE, corr = FALSE, method = "BFGS", warn.inf = TRUE) {

  nllbvpct.sym <- function(scale1, shape1, scale2, shape2, alpha) {
    nllbvpct(scale1, shape1, scale2, shape2, alpha, alpha)
  }

  nllbvpct <- function(scale1, shape1, scale2, shape2, alpha, beta) {
    .C("nllbvpct", spx$x1, spx$x2, spx$nn, spx$n, spx$thdi, spx$r1, spx$r2,
      spx$lambda, alpha, beta, scale1, shape1, scale2, shape2,
      dns = double(1), PACKAGE = "evd")$dns
  }
  param <- c("scale1", "shape1", "scale2", "shape2", "alpha")
  if(!sym) param <- c(param, "beta")
  nmdots <- names(list(...))
  start <- bvstart.vals(x, start, NULL, NULL, nmdots, param, NULL,
    NULL, model = "ct", sym = sym, obj = "bvpot", u = u)
  spx <- sep.bvdata(x, obj = "bvpot", u = u, censored = FALSE)
  nm <- names(start)
  fixed.param <- list(...)[nmdots %in% param]
  if(any(!(param %in% c(nm,names(fixed.param)))))
    stop("unspecified parameters")
  l <- length(nm)
  if(sym) f <- formals(nllbvpct.sym)
  else f <- formals(nllbvpct)
  names(f) <- param
  m <- match(nm, param)
  if(any(is.na(m))) 
    stop("`start' specifies unknown arguments") 
  if(sym) {
    formals(nllbvpct.sym) <- c(f[m], f[-m])
    nll <- function(p, ...) nllbvpct.sym(p, ...)
    if(l > 1) {
      body(nll) <- parse(text = paste("nllbvpct.sym(", paste("p[", 1:l, "]", collapse=", "), ", ...)"))
    }
  } else {
    formals(nllbvpct) <- c(f[m], f[-m])
    nll <- function(p, ...) nllbvpct(p, ...)
    if(l > 1) {
      body(nll) <- parse(text = paste("nllbvpct(", paste("p[", 1:l, "]", collapse=", "), ", ...)"))
    }
  }
  start.arg <- c(list(p = unlist(start)), fixed.param)
  if(warn.inf && do.call("nll", start.arg) == 1e+06)
    warning("negative log-likelihood is infinite at starting values")
  opt <- optim(start, nll, hessian = std.err, ..., method = method)
  if(sym) {
    beta <- opt$par["alpha"]
    names(beta) <- NULL
    if(is.na(beta)) beta <- fixed.param[["alpha"]]
    fixed.param <- c(fixed.param, list(beta = beta))
  }
  bvtpost.optim(x, u, opt, nm, fixed.param, std.err, dsm, corr, spx$nat, sym, model = "ct")
}

fbvpbilog <- function(x, u, start, ..., std.err = TRUE, dsm = TRUE, corr = FALSE, method = "BFGS", warn.inf = TRUE) {

  nllbvpbilog <- function(scale1, shape1, scale2, shape2, alpha, beta) {
    .C("nllbvpbilog", spx$x1, spx$x2, spx$nn, spx$n, spx$thdi, spx$r1,
      spx$r2, spx$lambda, alpha, beta, scale1, shape1, scale2, shape2,
      dns = double(1), PACKAGE = "evd")$dns    
  }
  param <- c("scale1", "shape1", "scale2", "shape2", "alpha", "beta")
  nmdots <- names(list(...))
  start <- bvstart.vals(x, start, NULL, NULL, nmdots, param, NULL,
    NULL, model = "bilog", sym = sym, obj = "bvpot", u = u)
  spx <- sep.bvdata(x, obj = "bvpot", u = u, censored = FALSE) 
  nm <- names(start)
  fixed.param <- list(...)[nmdots %in% param]
  if(any(!(param %in% c(nm,names(fixed.param)))))
    stop("unspecified parameters")
  l <- length(nm)
  f <- formals(nllbvpbilog)
  names(f) <- param
  m <- match(nm, param)
  if(any(is.na(m))) 
    stop("`start' specifies unknown arguments") 
  formals(nllbvpbilog) <- c(f[m], f[-m])
  nll <- function(p, ...) nllbvpbilog(p, ...)
  if(l > 1) {
    body(nll) <- parse(text = paste("nllbvpbilog(", paste("p[", 1:l,
      "]", collapse=", "), ", ...)"))
  }
  start.arg <- c(list(p = unlist(start)), fixed.param)
  if(warn.inf && do.call("nll", start.arg) == 1e+06)
    warning("negative log-likelihood is infinite at starting values")
  opt <- optim(start, nll, hessian = std.err, ..., method = method)
  bvtpost.optim(x, u, opt, nm, fixed.param, std.err, dsm, corr, spx$nat, FALSE, model = "bilog")
}

fbvpnegbilog <- function(x, u, start, ..., std.err = TRUE, dsm = TRUE, corr = FALSE, method = "BFGS", warn.inf = TRUE) {

  nllbvpnegbilog <- function(scale1, shape1, scale2, shape2, alpha, beta) {
    .C("nllbvpnegbilog", spx$x1, spx$x2, spx$nn, spx$n, spx$thdi, spx$r1,
      spx$r2, spx$lambda, alpha, beta, scale1, shape1, scale2, shape2,
      dns = double(1), PACKAGE = "evd")$dns 
  }
  param <- c("scale1", "shape1", "scale2", "shape2", "alpha", "beta")
  nmdots <- names(list(...))
  start <- bvstart.vals(x, start, NULL, NULL, nmdots, param, NULL,
    NULL, model = "negbilog", sym = sym, obj = "bvpot", u = u)
  spx <- sep.bvdata(x, obj = "bvpot", u = u, censored = FALSE) 
  nm <- names(start)
  fixed.param <- list(...)[nmdots %in% param]
  if(any(!(param %in% c(nm,names(fixed.param)))))
    stop("unspecified parameters")
  l <- length(nm)
  f <- formals(nllbvpnegbilog)
  names(f) <- param
  m <- match(nm, param)
  if(any(is.na(m))) 
    stop("`start' specifies unknown arguments") 
  formals(nllbvpnegbilog) <- c(f[m], f[-m])
  nll <- function(p, ...) nllbvpnegbilog(p, ...)
  if(l > 1) {
    body(nll) <- parse(text = paste("nllbvpnegbilog(", paste("p[", 1:l,
      "]", collapse=", "), ", ...)"))
  }
  start.arg <- c(list(p = unlist(start)), fixed.param)
  if(warn.inf && do.call("nll", start.arg) == 1e+06)
    warning("negative log-likelihood is infinite at starting values")
  opt <- optim(start, nll, hessian = std.err, ..., method = method)
  bvtpost.optim(x, u, opt, nm, fixed.param, std.err, dsm, corr, spx$nat, FALSE, model = "negbilog")
}


###################### post-optimisation processing #####################

bvtpost.optim <- function(x, u, opt, nm, fixed.param, std.err, dsm, corr, nat, sym, model) {
  if(opt$convergence != 0) {
    warning(paste("optimization for", model, "may not have succeeded"), call. = FALSE)
      if(opt$convergence == 1) 
        opt$convergence <- "iteration limit reached"
  }
  else opt$convergence <- "successful"
  if(std.err) {
    tol <- .Machine$double.eps^0.5
    var.cov <- qr(opt$hessian, tol = tol)
    if(var.cov$rank != ncol(var.cov$qr)) 
      stop(paste("observed information matrix for", model, "is singular; use std.err = FALSE"))
    var.cov <- solve(var.cov, tol = tol)
    std.err <- diag(var.cov)
    if(any(std.err <= 0)) 
      stop(paste("observed information matrix for", model, "is singular; use std.err = FALSE"))
    std.err <- sqrt(std.err)
    names(std.err) <- nm
    if(corr) {
      .mat <- diag(1/std.err, nrow = length(std.err))
      corr <- structure(.mat %*% var.cov %*% .mat, dimnames = list(nm, nm))
      diag(corr) <- rep(1, length(std.err))
    }
    else corr <- NULL
  }
  else std.err <- corr <- NULL
  param <- c(opt$par, unlist(fixed.param))
  if(dsm) {
    dep.sum <- numeric(3)
    if(model %in% c("log", "hr", "neglog")) {
      dep <- param["dep"]
      dep.sum[1] <- 2 * (1 - abvpar(dep = dep, model = model))
      dep.sum[2] <- 4 * integrate(function(x) 1 - abvpar(x, dep = dep, model = model), 0, 1)$value
    }
    if(model %in% c("alog", "aneglog")) {
      dep <- param["dep"]
      asy <- param[c("asy1", "asy2")]
      dep.sum[1] <- 2 * (1 - abvpar(dep = dep, asy = asy, model = model))
      dep.sum[2] <- 4 * integrate(function(x) 1 - abvpar(x, dep = dep, asy = asy, model = model), 0, 1)$value
      dffn <- function(x) abvpar(x, dep = dep, asy = asy, model = model) - abvpar(x, dep = dep, asy = rev(asy), model = model)
      dep.sum[3] <- 4 * integrate(dffn, 0, 0.5)$value/(3 - 2 * sqrt(2))
    }
    if(model %in% c("bilog", "negbilog", "ct")) {
      alpha <- param["alpha"]
      beta <- param["beta"]
      dep.sum[1] <- 2 * (1 - abvpar(alpha = alpha, beta = beta, model = model))
      dep.sum[2] <- 4 * integrate(function(x) 1 - abvpar(x, alpha = alpha, beta = beta, model = model), 0, 1)$value
      dffn <- function(x) abvpar(x, alpha = alpha, beta = beta, model = model) - abvpar(x, alpha = beta, beta = alpha, model = model)
      dep.sum[3] <- 4 * integrate(dffn, 0, 0.5)$value/(3 - 2 * sqrt(2))
    }
  }
  else dep.sum <- NULL
  list(estimate = opt$par, std.err = std.err, fixed = unlist(fixed.param), param = param, deviance = 2 * opt$value, dep.summary = dep.sum, corr = corr, convergence = opt$convergence, counts = opt$counts, message = opt$message, data = x, threshold = u, n = nrow(x), nat = nat, sym = sym, model = model)
}


########################## method functions #########################


"print.bvpot" <-  function(x, digits = max(3, getOption("digits") - 3), ...) 
{
    cat("\nCall:", deparse(x$call), "\n")
    cat("Deviance:", deviance(x), "\n")
    cat("AIC:", AIC(x), "\n")

    cat("\nThreshold:", round(x$threshold, digits), "\n")
    cat("Marginal Number Above:", x$nat[1:2], "\n")
    cat("Marginal Proportion Above:", round(x$nat[1:2]/x$n, digits), "\n")
    cat("Number Above:", x$nat[3], "\n")
    cat("Proportion Above:", round(x$nat[3]/x$n, digits), "\n")
    
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

plot.bvpot <- function(x, ...) {
  stop("no plot method currently implemented for bivariate threshold models")
}
