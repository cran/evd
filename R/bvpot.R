
################## bivariate threshold fitting routines #################

fbvpot <- function(x, threshold, model = c("log", "bilog", "alog", "neglog", "negbilog", "aneglog", "ct", "hr"), likelihood = c("censored","poisson"), start, ..., sym = FALSE, cshape = cscale, cscale = FALSE, std.err = TRUE, dsm = TRUE, corr = FALSE, method = "BFGS", warn.inf = TRUE) {
  call <- match.call()
  likelihood <- match.arg(likelihood)
  ft <- switch(likelihood,
    censored = fbvcpot(x = x, u = threshold, model = model, start = start, ...,
      sym = sym, cshape = cshape, cscale = cscale, std.err =
      std.err, dsm = dsm, corr = corr, method = method, warn.inf = warn.inf),
    poisson = fbvppot(x = x, u = threshold, model = model, start = start, ...,
      sym = sym, cshape = cshape, cscale = cscale, std.err =
      std.err, dsm = dsm, corr = corr, method = method, warn.inf = warn.inf))
  structure(c(ft, call = call), class = c("bvpot", "evd"))
}

fbvcpot <- function(x, u, model = c("log", "bilog", "alog", "neglog", "negbilog", "aneglog", "ct", "hr"), start, ..., sym = FALSE, cshape = cscale, cscale = FALSE, std.err = TRUE, dsm = TRUE, corr = FALSE, method = "BFGS", warn.inf = TRUE) {
  model <- match.arg(model)
  if(sym && !(model %in% c("alog","aneglog","ct")))
    warning("Argument `sym' was ignored")
  switch(model,
    log = fbvclog(x = x, u = u, start = start, ..., cshape = cshape,
      cscale = cscale, std.err = std.err, dsm = dsm,
      corr = corr, method = method, warn.inf = warn.inf),
    bilog = fbvcbilog(x = x, u = u, start = start, ..., cshape = cshape,
      cscale = cscale, std.err = std.err, dsm = dsm,
      corr = corr, method = method, warn.inf = warn.inf),
    alog = fbvcalog(x = x, u = u, start = start, ..., sym = sym,
      cshape = cshape, cscale = cscale, std.err = std.err,
      dsm = dsm, corr = corr, method = method, warn.inf = warn.inf),
    neglog = fbvcneglog(x = x, u = u, start = start, ..., cshape =
      cshape, cscale = cscale, std.err = std.err,
      dsm = dsm, corr = corr, method = method, warn.inf = warn.inf),
    negbilog = fbvcnegbilog(x = x, u = u, start = start, ..., cshape =
      cshape, cscale = cscale, std.err = std.err, dsm = dsm,
      corr = corr, method = method, warn.inf = warn.inf),
    aneglog = fbvcaneglog(x = x, u = u, start = start, ..., sym = sym,
      cshape = cshape, cscale = cscale, std.err = std.err,
      dsm = dsm, corr = corr, method = method, warn.inf = warn.inf),
    ct = fbvcct(x = x, u = u, start = start, ..., sym = sym, cshape =
      cshape, cscale = cscale, std.err = std.err, dsm = dsm,
      corr = corr, method = method, warn.inf = warn.inf),
    hr = fbvchr(x = x, u = u, start = start, ..., cshape = cshape,
      cscale = cscale, std.err = std.err, dsm = dsm, corr =
      corr, method = method, warn.inf = warn.inf))
}

fbvppot <- function(x, u, model = c("log", "bilog", "alog", "neglog", "negbilog", "aneglog", "ct", "hr"), start, ..., sym = FALSE, cshape = cscale, cscale = FALSE, std.err = TRUE, dsm = TRUE, corr = FALSE, method = "BFGS", warn.inf = TRUE) {
  model <- match.arg(model)
  if(model %in% c("alog","aneglog","hr"))
    stop("This model is not implemented for poisson likelihood")
  if(sym && (model != "ct"))
    warning("Argument `sym' was ignored")
  switch(model,
    log = fbvplog(x = x, u = u, start = start, ..., cshape = cshape,
      cscale = cscale, std.err = std.err, dsm = dsm,
      corr = corr, method = method, warn.inf = warn.inf),
    bilog = fbvpbilog(x = x, u = u, start = start, ..., cshape = cshape,
      cscale = cscale, std.err = std.err, dsm = dsm, corr =
      corr, method = method, warn.inf = warn.inf),
    neglog = fbvpneglog(x = x, u = u, start = start, ..., cshape =
      cshape, cscale = cscale, std.err = std.err,
      dsm = dsm, corr = corr, method = method, warn.inf = warn.inf),
    negbilog = fbvpnegbilog(x = x, u = u, start = start, ..., cshape =
      cshape, cscale = cscale, std.err = std.err, dsm = dsm,
      corr = corr, method = method, warn.inf = warn.inf),
    ct = fbvpct(x = x, u = u, start = start, ..., sym = sym, cshape =
      cshape, cscale = cscale, std.err = std.err, dsm = dsm,
      corr = corr, method = method, warn.inf = warn.inf))
}

################## censored likelihood fitting routines #################

fbvclog <- function(x, u, start, ..., cshape = cscale, cscale = FALSE, std.err = TRUE, dsm = TRUE, corr = FALSE, method = "BFGS", warn.inf = TRUE) {

  nllbvclog <- function(scale1, shape1, scale2, shape2, dep) {
    if(cshape) shape2 <- shape1
    if(cscale) scale2 <- scale1
    .C("nllbvclog", spx$x1, spx$x2, spx$nn, spx$n, spx$thdi,
      spx$lambda, dep, scale1, shape1, scale2, shape2, dns = double(1),
      PACKAGE = "evd")$dns
  }
  param <- c("scale1", "shape1")
  if(!cscale) param <- c(param, "scale2")
  if(!cshape) param <- c(param, "shape2")
  param <- c(param, "dep")
  nmdots <- names(list(...))
  start <- bvstart.vals(x, start, NULL, NULL, nmdots, param, NULL,
    NULL, model = "log", obj = "bvpot", u = u)
  spx <- sep.bvdata(x, obj = "bvpot", u = u)  
  nm <- names(start)
  l <- length(nm)
  fixed.param <- list(...)[nmdots %in% param]
  if(any(!(param %in% c(nm,names(fixed.param)))))
    stop("unspecified parameters")
  prind <- c(TRUE, TRUE, !cscale, !cshape, TRUE)
  f <- formals(nllbvclog)[prind]
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
  bvtpost.optim(x, u, opt, nm, fixed.param, std.err, dsm, corr, spx$nat, sym = FALSE, cmar = c(cscale, cshape), model = "log")
}

fbvcbilog <- function(x, u, start, ..., cshape = cscale, cscale = FALSE, std.err = TRUE, dsm = TRUE, corr = FALSE, method = "BFGS", warn.inf = TRUE) {

  nllbvcbilog <- function(scale1, shape1, scale2, shape2, alpha, beta) {
    if(cshape) shape2 <- shape1
    if(cscale) scale2 <- scale1
    .C("nllbvcbilog", spx$x1, spx$x2, spx$nn, spx$n, spx$thdi, spx$lambda,
      alpha, beta, scale1, shape1, scale2, shape2, dns = double(1),
      PACKAGE = "evd")$dns
  }
  param <- c("scale1", "shape1")
  if(!cscale) param <- c(param, "scale2")
  if(!cshape) param <- c(param, "shape2")
  param <- c(param, "alpha", "beta")
  nmdots <- names(list(...))
  start <- bvstart.vals(x, start, NULL, NULL, nmdots, param, NULL,
    NULL, model = "bilog", obj = "bvpot", u = u)
  spx <- sep.bvdata(x, obj = "bvpot", u = u) 
  nm <- names(start)
  l <- length(nm)
  fixed.param <- list(...)[nmdots %in% param]
  if(any(!(param %in% c(nm,names(fixed.param)))))
    stop("unspecified parameters")
  prind <- c(TRUE, TRUE, !cscale, !cshape, TRUE, TRUE)
  f <- formals(nllbvcbilog)[prind]
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
  bvtpost.optim(x, u, opt, nm, fixed.param, std.err, dsm, corr, spx$nat, sym = FALSE, cmar = c(cscale, cshape), model = "bilog")
}

fbvcalog <- function(x, u, start, ..., sym = FALSE, cshape = cscale, cscale = FALSE, std.err = TRUE, dsm = TRUE, corr = FALSE, method = "BFGS", warn.inf = TRUE) {
  nllbvcalog <- function(scale1, shape1, scale2, shape2, asy1, asy2, dep) {
    if(sym) asy2 <- asy1
    if(cshape) shape2 <- shape1
    if(cscale) scale2 <- scale1
    .C("nllbvcalog", spx$x1, spx$x2, spx$nn, spx$n, spx$thdi, spx$lambda,
       dep, asy1, asy2, scale1, shape1, scale2, shape2, dns = double(1),
       PACKAGE = "evd")$dns
  }
  param <- c("scale1", "shape1")
  if(!cscale) param <- c(param, "scale2")
  if(!cshape) param <- c(param, "shape2")
  if(!sym) param <- c(param, "asy1", "asy2", "dep")
  else param <- c(param, "asy1", "dep")
  nmdots <- names(list(...))
  start <- bvstart.vals(x, start, NULL, NULL, nmdots, param, NULL,
    NULL, model = "alog", obj = "bvpot", u = u)
  spx <- sep.bvdata(x, obj = "bvpot", u = u)
  nm <- names(start)
  l <- length(nm)
  fixed.param <- list(...)[nmdots %in% param]
  if(any(!(param %in% c(nm,names(fixed.param)))))
    stop("unspecified parameters")
  prind <- c(TRUE, TRUE, !cscale, !cshape, TRUE, !sym, TRUE)
  f <- formals(nllbvcalog)[prind]
  names(f) <- param
  m <- match(nm, param)
  if(any(is.na(m))) 
    stop("`start' specifies unknown arguments") 
  formals(nllbvcalog) <- c(f[m], f[-m])
  nll <- function(p, ...) nllbvcalog(p, ...)
  if(l > 1) {
    body(nll) <- parse(text = paste("nllbvcalog(", paste("p[", 1:l, "]", collapse=", "), ", ...)"))
  }
  start.arg <- c(list(p = unlist(start)), fixed.param)
  if(warn.inf && do.call("nll", start.arg) == 1e+06)
    warning("negative log-likelihood is infinite at starting values")
  opt <- optim(start, nll, hessian = std.err, ..., method = method)
  bvtpost.optim(x, u, opt, nm, fixed.param, std.err, dsm, corr,
    spx$nat, sym = sym, cmar = c(cscale, cshape), model = "alog")
}

fbvcneglog <- function(x, u, start, ..., cshape = cscale, cscale = FALSE, std.err = TRUE, dsm = TRUE, corr = FALSE, method = "BFGS", warn.inf = TRUE) {

  nllbvcneglog <- function(scale1, shape1, scale2, shape2, dep) {
    if(cshape) shape2 <- shape1
    if(cscale) scale2 <- scale1
    .C("nllbvcneglog", spx$x1, spx$x2, spx$nn, spx$n, spx$thdi, spx$lambda,
      dep, scale1, shape1, scale2, shape2, dns = double(1),
      PACKAGE = "evd")$dns
  }
  param <- c("scale1", "shape1")
  if(!cscale) param <- c(param, "scale2")
  if(!cshape) param <- c(param, "shape2")
  param <- c(param, "dep")
  nmdots <- names(list(...))
  start <- bvstart.vals(x, start, NULL, NULL, nmdots, param, NULL,
    NULL, model = "neglog", obj = "bvpot", u = u)
  spx <- sep.bvdata(x, obj = "bvpot", u = u) 
  nm <- names(start)
  l <- length(nm)
  fixed.param <- list(...)[nmdots %in% param]
  if(any(!(param %in% c(nm,names(fixed.param)))))
    stop("unspecified parameters")
  prind <- c(TRUE, TRUE, !cscale, !cshape, TRUE)
  f <- formals(nllbvcneglog)[prind]
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
  bvtpost.optim(x, u, opt, nm, fixed.param, std.err, dsm, corr, spx$nat, sym = FALSE, cmar = c(cscale, cshape), model = "neglog")
}

fbvcnegbilog <- function(x, u, start, ..., cshape = cscale, cscale = FALSE, std.err = TRUE, dsm = TRUE, corr = FALSE, method = "BFGS", warn.inf = TRUE) {

  nllbvcnegbilog <- function(scale1, shape1, scale2, shape2, alpha, beta) {
    if(cshape) shape2 <- shape1
    if(cscale) scale2 <- scale1
    .C("nllbvcnegbilog", spx$x1, spx$x2, spx$nn, spx$n, spx$thdi, spx$lambda,
       alpha, beta, scale1, shape1, scale2, shape2, dns = double(1),
       PACKAGE = "evd")$dns
  }
  param <- c("scale1", "shape1")
  if(!cscale) param <- c(param, "scale2")
  if(!cshape) param <- c(param, "shape2")
  param <- c(param, "alpha", "beta")
  nmdots <- names(list(...))
  start <- bvstart.vals(x, start, NULL, NULL, nmdots, param, NULL,
    NULL, model = "negbilog", obj = "bvpot", u = u)
  spx <- sep.bvdata(x, obj = "bvpot", u = u) 
  nm <- names(start)
  l <- length(nm)
  fixed.param <- list(...)[nmdots %in% param]
  if(any(!(param %in% c(nm,names(fixed.param)))))
    stop("unspecified parameters")
  prind <- c(TRUE, TRUE, !cscale, !cshape, TRUE, TRUE)
  f <- formals(nllbvcnegbilog)[prind]
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
  bvtpost.optim(x, u, opt, nm, fixed.param, std.err, dsm, corr, spx$nat, sym = FALSE, cmar = c(cscale, cshape), model = "negbilog")
}

fbvcaneglog <- function(x, u, start, ..., sym = FALSE, cshape = cscale, cscale = FALSE, std.err = TRUE, dsm = TRUE, corr = FALSE, method = "BFGS", warn.inf = TRUE) {
  
  nllbvcaneglog <- function(scale1, shape1, scale2, shape2, asy1, asy2, dep) {
    if(sym) asy2 <- asy1
    if(cshape) shape2 <- shape1
    if(cscale) scale2 <- scale1
    .C("nllbvcaneglog", spx$x1, spx$x2, spx$nn, spx$n, spx$thdi, spx$lambda,
       dep, asy1, asy2, scale1, shape1, scale2, shape2, dns = double(1),
       PACKAGE = "evd")$dns
  }
  param <- c("scale1", "shape1")
  if(!cscale) param <- c(param, "scale2")
  if(!cshape) param <- c(param, "shape2")
  if(!sym) param <- c(param, "asy1", "asy2", "dep")
  else param <- c(param, "asy1", "dep")
  nmdots <- names(list(...))
  start <- bvstart.vals(x, start, NULL, NULL, nmdots, param, NULL,
    NULL, model = "aneglog", obj = "bvpot", u = u)
  spx <- sep.bvdata(x, obj = "bvpot", u = u) 
  nm <- names(start)
  l <- length(nm)
  fixed.param <- list(...)[nmdots %in% param]
  if(any(!(param %in% c(nm,names(fixed.param)))))
    stop("unspecified parameters")
  prind <- c(TRUE, TRUE, !cscale, !cshape, TRUE, !sym, TRUE)
  f <- formals(nllbvcaneglog)[prind]
  names(f) <- param
  m <- match(nm, param)
  if(any(is.na(m))) 
    stop("`start' specifies unknown arguments") 
  formals(nllbvcaneglog) <- c(f[m], f[-m])
  nll <- function(p, ...) nllbvcaneglog(p, ...)
  if(l > 1) {
    body(nll) <- parse(text = paste("nllbvcaneglog(", paste("p[", 1:l, "]", collapse=", "), ", ...)"))
  }
  start.arg <- c(list(p = unlist(start)), fixed.param)
  if(warn.inf && do.call("nll", start.arg) == 1e+06)
    warning("negative log-likelihood is infinite at starting values")
  opt <- optim(start, nll, hessian = std.err, ..., method = method)
  bvtpost.optim(x, u, opt, nm, fixed.param, std.err, dsm, corr, spx$nat, sym = sym, cmar = c(cscale, cshape), model = "aneglog")
}

fbvcct <- function(x, u, start, ..., sym = FALSE, cshape = cscale, cscale = FALSE, std.err = TRUE, dsm = TRUE, corr = FALSE, method = "BFGS", warn.inf = TRUE) {

  nllbvcct <- function(scale1, shape1, scale2, shape2, alpha, beta) {
    if(sym) beta <- alpha
    if(cshape) shape2 <- shape1
    if(cscale) scale2 <- scale1
    .C("nllbvcct", spx$x1, spx$x2, spx$nn, spx$n, spx$thdi, spx$lambda,
      alpha, beta, scale1, shape1, scale2, shape2, dns = double(1),
      PACKAGE = "evd")$dns
  }
  param <- c("scale1", "shape1")
  if(!cscale) param <- c(param, "scale2")
  if(!cshape) param <- c(param, "shape2")
  if(!sym) param <- c(param, "alpha", "beta")
  else param <- c(param, "alpha")
  nmdots <- names(list(...))
  start <- bvstart.vals(x, start, NULL, NULL, nmdots, param, NULL,
    NULL, model = "ct", obj = "bvpot", u = u)
  spx <- sep.bvdata(x, obj = "bvpot", u = u) 
  nm <- names(start)
  l <- length(nm)
  fixed.param <- list(...)[nmdots %in% param]
  if(any(!(param %in% c(nm,names(fixed.param)))))
    stop("unspecified parameters")
  prind <- c(TRUE, TRUE, !cscale, !cshape, TRUE, !sym)
  f <- formals(nllbvcct)[prind]
  names(f) <- param
  m <- match(nm, param)
  if(any(is.na(m))) 
    stop("`start' specifies unknown arguments") 
  formals(nllbvcct) <- c(f[m], f[-m])
  nll <- function(p, ...) nllbvcct(p, ...)
  if(l > 1) {
    body(nll) <- parse(text = paste("nllbvcct(", paste("p[", 1:l, "]", collapse=", "), ", ...)"))
  }
  start.arg <- c(list(p = unlist(start)), fixed.param)
  if(warn.inf && do.call("nll", start.arg) == 1e+06)
    warning("negative log-likelihood is infinite at starting values")
  opt <- optim(start, nll, hessian = std.err, ..., method = method)
  bvtpost.optim(x, u, opt, nm, fixed.param, std.err, dsm, corr, spx$nat, sym = sym, cmar = c(cscale, cshape), model = "ct")
}

fbvchr <- function(x, u, start, ..., cshape = cscale, cscale = FALSE, std.err = TRUE, dsm = TRUE, corr = FALSE, method = "BFGS", warn.inf = TRUE) {

  nllbvchr <- function(scale1, shape1, scale2, shape2, dep) {
    if(cshape) shape2 <- shape1
    if(cscale) scale2 <- scale1
    .C("nllbvchr", spx$x1, spx$x2, spx$nn, spx$n, spx$thdi,
      spx$lambda, dep, scale1, shape1, scale2, shape2, dns = double(1),
      PACKAGE = "evd")$dns
  }
  param <- c("scale1", "shape1")
  if(!cscale) param <- c(param, "scale2")
  if(!cshape) param <- c(param, "shape2")
  param <- c(param, "dep")
  nmdots <- names(list(...))
  start <- bvstart.vals(x, start, NULL, NULL, nmdots, param, NULL,
    NULL, model = "hr", obj = "bvpot", u = u)
  spx <- sep.bvdata(x, obj = "bvpot", u = u)  
  nm <- names(start)
  l <- length(nm)
  fixed.param <- list(...)[nmdots %in% param]
  if(any(!(param %in% c(nm,names(fixed.param)))))
    stop("unspecified parameters")
  prind <- c(TRUE, TRUE, !cscale, !cshape, TRUE)
  f <- formals(nllbvchr)[prind]
  names(f) <- param
  m <- match(nm, param)
  if(any(is.na(m))) 
    stop("`start' specifies unknown arguments")     
  formals(nllbvchr) <- c(f[m], f[-m])
  nll <- function(p, ...) nllbvchr(p, ...)
  if(l > 1) {
    body(nll) <- parse(text = paste("nllbvchr(", paste("p[", 1:l, "]", collapse=", "), ", ...)"))
  }
  start.arg <- c(list(p = unlist(start)), fixed.param)
  if(warn.inf && do.call("nll", start.arg) == 1e+06)
    warning("negative log-likelihood is infinite at starting values")
  opt <- optim(start, nll, hessian = std.err, ..., method = method)
  bvtpost.optim(x, u, opt, nm, fixed.param, std.err, dsm, corr, spx$nat, sym = FALSE, cmar = c(cscale, cshape), model = "hr")
}

################## Poisson likelihood fitting routines ##################

fbvplog <- function(x, u, start, ..., cshape = cscale, cscale = FALSE, std.err = TRUE, dsm = TRUE, corr = FALSE, method = "BFGS", warn.inf = TRUE) {

  nllbvplog <- function(scale1, shape1, scale2, shape2, dep) {
    if(cshape) shape2 <- shape1
    if(cscale) scale2 <- scale1
    .C("nllbvplog", spx$x1, spx$x2, spx$nn, spx$n, spx$thdi, spx$r1, spx$r2,
      spx$lambda, dep, scale1, shape1, scale2, shape2, dns = double(1),
      PACKAGE = "evd")$dns
  }
  param <- c("scale1", "shape1")
  if(!cscale) param <- c(param, "scale2")
  if(!cshape) param <- c(param, "shape2")
  param <- c(param, "dep")
  nmdots <- names(list(...))
  start <- bvstart.vals(x, start, NULL, NULL, nmdots, param, NULL,
    NULL, model = "log", obj = "bvpot", u = u)
  spx <- sep.bvdata(x, obj = "bvpot", u = u, censored = FALSE) 
  nm <- names(start)
  l <- length(nm)
  fixed.param <- list(...)[nmdots %in% param]
  if(any(!(param %in% c(nm,names(fixed.param)))))
    stop("unspecified parameters")
  prind <- c(TRUE, TRUE, !cscale, !cshape, TRUE)
  f <- formals(nllbvplog)[prind]
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
    spx$nat, sym = FALSE, cmar = c(cscale, cshape), model = "log")
}

fbvpneglog <- function(x, u, start, ..., cshape = cscale, cscale = FALSE, std.err = TRUE, dsm = TRUE, corr = FALSE, method = "BFGS", warn.inf = TRUE) {

  nllbvpneglog <- function(scale1, shape1, scale2, shape2, dep) {
    if(cshape) shape2 <- shape1
    if(cscale) scale2 <- scale1
    .C("nllbvpneglog", spx$x1, spx$x2, spx$nn, spx$n, spx$thdi, spx$r1,
      spx$r2, spx$lambda, dep, scale1, shape1, scale2, shape2,
      dns = double(1), PACKAGE = "evd")$dns    
  }
  param <- c("scale1", "shape1")
  if(!cscale) param <- c(param, "scale2")
  if(!cshape) param <- c(param, "shape2")
  param <- c(param, "dep")
  nmdots <- names(list(...))
  start <- bvstart.vals(x, start, NULL, NULL, nmdots, param, NULL,
    NULL, model = "neglog", obj = "bvpot", u = u)
  spx <- sep.bvdata(x, obj = "bvpot", u = u, censored = FALSE)
  nm <- names(start)
  l <- length(nm)
  fixed.param <- list(...)[nmdots %in% param]
  if(any(!(param %in% c(nm,names(fixed.param)))))
    stop("unspecified parameters")
  prind <- c(TRUE, TRUE, !cscale, !cshape, TRUE)
  f <- formals(nllbvpneglog)[prind]
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
  bvtpost.optim(x, u, opt, nm, fixed.param, std.err, dsm, corr, spx$nat, sym = FALSE, cmar = c(cscale, cshape), model = "neglog")
}

fbvpct <- function(x, u, start, ..., sym = FALSE, cshape = cscale, cscale = FALSE, std.err = TRUE, dsm = TRUE, corr = FALSE, method = "BFGS", warn.inf = TRUE) {

  nllbvpct <- function(scale1, shape1, scale2, shape2, alpha, beta) {
    if(sym) beta <- alpha
    if(cshape) shape2 <- shape1
    if(cscale) scale2 <- scale1
    .C("nllbvpct", spx$x1, spx$x2, spx$nn, spx$n, spx$thdi, spx$r1, spx$r2,
      spx$lambda, alpha, beta, scale1, shape1, scale2, shape2,
      dns = double(1), PACKAGE = "evd")$dns
  }
  param <- c("scale1", "shape1")
  if(!cscale) param <- c(param, "scale2")
  if(!cshape) param <- c(param, "shape2")
  if(!sym) param <- c(param, "alpha", "beta")
  else param <- c(param, "alpha")
  nmdots <- names(list(...))
  start <- bvstart.vals(x, start, NULL, NULL, nmdots, param, NULL,
    NULL, model = "ct", obj = "bvpot", u = u)
  spx <- sep.bvdata(x, obj = "bvpot", u = u, censored = FALSE)
  nm <- names(start)
  l <- length(nm)
  fixed.param <- list(...)[nmdots %in% param]
  if(any(!(param %in% c(nm,names(fixed.param)))))
    stop("unspecified parameters")
  prind <- c(TRUE, TRUE, !cscale, !cshape, TRUE, !sym)
  f <- formals(nllbvpct)[prind]
  names(f) <- param
  m <- match(nm, param)
  if(any(is.na(m))) 
    stop("`start' specifies unknown arguments") 
  formals(nllbvpct) <- c(f[m], f[-m])
  nll <- function(p, ...) nllbvpct(p, ...)
  if(l > 1) {
    body(nll) <- parse(text = paste("nllbvpct(", paste("p[", 1:l, "]", collapse=", "), ", ...)"))
  }
  start.arg <- c(list(p = unlist(start)), fixed.param)
  if(warn.inf && do.call("nll", start.arg) == 1e+06)
    warning("negative log-likelihood is infinite at starting values")
  opt <- optim(start, nll, hessian = std.err, ..., method = method)
  bvtpost.optim(x, u, opt, nm, fixed.param, std.err, dsm, corr, spx$nat, sym = sym, cmar = c(cscale, cshape), model = "ct")
}

fbvpbilog <- function(x, u, start, ..., cshape = cscale, cscale = FALSE, std.err = TRUE, dsm = TRUE, corr = FALSE, method = "BFGS", warn.inf = TRUE) {

  nllbvpbilog <- function(scale1, shape1, scale2, shape2, alpha, beta) {
    if(cshape) shape2 <- shape1
    if(cscale) scale2 <- scale1
    .C("nllbvpbilog", spx$x1, spx$x2, spx$nn, spx$n, spx$thdi, spx$r1,
      spx$r2, spx$lambda, alpha, beta, scale1, shape1, scale2, shape2,
      dns = double(1), PACKAGE = "evd")$dns    
  }
  param <- c("scale1", "shape1")
  if(!cscale) param <- c(param, "scale2")
  if(!cshape) param <- c(param, "shape2")
  param <- c(param, "alpha", "beta")
  nmdots <- names(list(...))
  start <- bvstart.vals(x, start, NULL, NULL, nmdots, param, NULL,
    NULL, model = "bilog", obj = "bvpot", u = u)
  spx <- sep.bvdata(x, obj = "bvpot", u = u, censored = FALSE) 
  nm <- names(start)
  l <- length(nm)
  fixed.param <- list(...)[nmdots %in% param]
  if(any(!(param %in% c(nm,names(fixed.param)))))
    stop("unspecified parameters")
  prind <- c(TRUE, TRUE, !cscale, !cshape, TRUE, TRUE)
  f <- formals(nllbvpbilog)[prind]
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
  bvtpost.optim(x, u, opt, nm, fixed.param, std.err, dsm, corr, spx$nat, sym = FALSE, cmar = c(cscale, cshape), model = "bilog")
}

fbvpnegbilog <- function(x, u, start, ..., cshape = cscale, cscale = FALSE, std.err = TRUE, dsm = TRUE, corr = FALSE, method = "BFGS", warn.inf = TRUE) {

  nllbvpnegbilog <- function(scale1, shape1, scale2, shape2, alpha, beta) {
    if(cshape) shape2 <- shape1
    if(cscale) scale2 <- scale1
    .C("nllbvpnegbilog", spx$x1, spx$x2, spx$nn, spx$n, spx$thdi, spx$r1,
      spx$r2, spx$lambda, alpha, beta, scale1, shape1, scale2, shape2,
      dns = double(1), PACKAGE = "evd")$dns 
  }
  param <- c("scale1", "shape1")
  if(!cscale) param <- c(param, "scale2")
  if(!cshape) param <- c(param, "shape2")
  param <- c(param, "alpha", "beta")
  nmdots <- names(list(...))
  start <- bvstart.vals(x, start, NULL, NULL, nmdots, param, NULL,
    NULL, model = "negbilog", obj = "bvpot", u = u)
  spx <- sep.bvdata(x, obj = "bvpot", u = u, censored = FALSE) 
  nm <- names(start)
  l <- length(nm)
  fixed.param <- list(...)[nmdots %in% param]
  if(any(!(param %in% c(nm,names(fixed.param)))))
    stop("unspecified parameters")
  prind <- c(TRUE, TRUE, !cscale, !cshape, TRUE, TRUE)
  f <- formals(nllbvpnegbilog)[prind]
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
  bvtpost.optim(x, u, opt, nm, fixed.param, std.err, dsm, corr, spx$nat, sym = FALSE, cmar = c(cscale, cshape), model = "negbilog")
}


###################### post-optimisation processing #####################

bvtpost.optim <- function(x, u, opt, nm, fixed.param, std.err, dsm, corr, nat, sym, cmar, model) {
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
  fixed <- unlist(fixed.param)
  param <- c(opt$par, fixed)
  fixed2 <- NULL
  if(cmar[1]) fixed2 <- c(fixed2, param["scale1"])
  if(cmar[2]) fixed2 <- c(fixed2, param["shape1"])
  if(sym) {
    if(model %in% c("alog","aneglog")) fixed2 <- c(fixed2, param["asy1"])
    if(model == "ct") fixed2 <- c(fixed2, param["alpha"])
  }
  if(!is.null(fixed2)) {
    names(fixed2) <- sub("1", "2", names(fixed2))
    names(fixed2) <- sub("alpha", "beta", names(fixed2))
  }
  param <- c(param, fixed2)
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
  list(estimate = opt$par, std.err = std.err, fixed = fixed, fixed2 = fixed2, param = param, deviance = 2 * opt$value, dep.summary = dep.sum, corr = corr, convergence = opt$convergence, counts = opt$counts, message = opt$message, data = x, threshold = u, n = nrow(x), nat = nat, sym = sym, cmar = cmar, model = model)
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

chiplot <- function(data, nq = 100, qlim = NULL, which = 1:2, conf = 0.95, lty = 1, cilty = 2, col = 1, cicol = 1, xlim = c(0,1), ylim1 = NULL, ylim2 = c(-1,1), main1 = "Chi Plot", main2 = "Chi Bar Plot", xlab = "Quantile", ylab1 = "Chi", ylab2 = "Chi Bar", ask = nb.fig < length(which) && dev.interactive(), ...)
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
    cnst <- qnorm((1 + conf)/2)
    varchi <- ((1/log(u)^2 * 1)/cu^2 * cu * (1 - cu))/n
    varchi <- cnst * sqrt(varchi)
    varchibar <- (((4 * log(1 - u)^2)/(log(cbaru)^4 * cbaru^2)) * cbaru * (
		1 - cbaru))/n
    varchibar <- cnst * sqrt(varchibar)

    chiu <- cbind(chilow = chiu-varchi, chi = chiu, chiupp = chiu+varchi) 
    chibaru <- cbind(chiblow = chibaru-varchibar, chib = chibaru, chibupp =
      chibaru+varchibar)

    show <- logical(2)
    show[which] <- TRUE
    lty <- c(cilty, lty, cilty)
    col <- c(cicol, col, cicol)
    nb.fig <- prod(par("mfcol"))
    if (ask) {
      op <- par(ask = TRUE)
      on.exit(par(op))
    }
    tmp <- rep(-Inf, nq)
    tmp[u > 0.5] <- 2 - log(2*u[u > 0.5] - 1)/log(u[u > 0.5])
    chiu <- pmin(pmax(chiu, tmp), 1)
    chibaru <- pmin(pmax(chibaru, -1), 1)
    if(is.null(ylim1)) ylim1 <- c(min(c(chiu, 0)), 1)
    if(show[1]) {
      matplot(u, chiu, type = "l", lty = lty, col = col, xlim = xlim, ylim = ylim1,
              main = main1, xlab = xlab, ylab = ylab1, ...)
    }

    if(show[2]) {
      matplot(u, chibaru, type = "l", lty = lty, col = col, xlim = xlim, ylim = ylim2,
              main = main2, xlab = xlab, ylab = ylab2, ...)
    }

    plvals <- list(quantile = u, chi = chiu, chibar = chibaru)
    if(!show[1]) plvals$chi <- NULL
    if(!show[2]) plvals$chib <- NULL
    invisible(plvals)    
}


