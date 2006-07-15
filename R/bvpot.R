################## bivariate threshold fitting routines #################

fbvpot <- function(x, threshold, model = c("log", "bilog", "alog", "neglog", "negbilog", "aneglog", "ct", "hr", "amix"), likelihood = c("censored","poisson"), start, ..., sym = FALSE, cshape = cscale, cscale = FALSE, std.err = TRUE, corr = FALSE, method = "BFGS", warn.inf = TRUE) {
  call <- match.call()
  likelihood <- match.arg(likelihood)
  ft <- switch(likelihood,
    censored = fbvcpot(x = x, u = threshold, model = model, start = start, ...,
      sym = sym, cshape = cshape, cscale = cscale, std.err =
      std.err, corr = corr, method = method, warn.inf = warn.inf),
    poisson = fbvppot(x = x, u = threshold, model = model, start = start, ...,
      sym = sym, cshape = cshape, cscale = cscale, std.err =
      std.err, corr = corr, method = method, warn.inf = warn.inf))
  structure(c(ft, call = call), class = c("bvpot", "evd"))
}

fbvcpot <- function(x, u, model = c("log", "bilog", "alog", "neglog", "negbilog", "aneglog", "ct", "hr", "amix"), start, ..., sym = FALSE, cshape = cscale, cscale = FALSE, std.err = TRUE, corr = FALSE, method = "BFGS", warn.inf = TRUE) {
  model <- match.arg(model)
  if(sym && !(model %in% c("alog","aneglog","ct")))
    warning("Argument `sym' was ignored")
  switch(model,
    log = fbvclog(x = x, u = u, start = start, ..., cshape = cshape,
      cscale = cscale, std.err = std.err,
      corr = corr, method = method, warn.inf = warn.inf),
    bilog = fbvcbilog(x = x, u = u, start = start, ..., cshape = cshape,
      cscale = cscale, std.err = std.err,
      corr = corr, method = method, warn.inf = warn.inf),
    alog = fbvcalog(x = x, u = u, start = start, ..., sym = sym,
      cshape = cshape, cscale = cscale, std.err = std.err,
      corr = corr, method = method, warn.inf = warn.inf),
    neglog = fbvcneglog(x = x, u = u, start = start, ..., cshape =
      cshape, cscale = cscale, std.err = std.err,
      corr = corr, method = method, warn.inf = warn.inf),
    negbilog = fbvcnegbilog(x = x, u = u, start = start, ..., cshape =
      cshape, cscale = cscale, std.err = std.err,
      corr = corr, method = method, warn.inf = warn.inf),
    aneglog = fbvcaneglog(x = x, u = u, start = start, ..., sym = sym,
      cshape = cshape, cscale = cscale, std.err = std.err,
      corr = corr, method = method, warn.inf = warn.inf),
    ct = fbvcct(x = x, u = u, start = start, ..., sym = sym, cshape =
      cshape, cscale = cscale, std.err = std.err,
      corr = corr, method = method, warn.inf = warn.inf),
    hr = fbvchr(x = x, u = u, start = start, ..., cshape = cshape,
      cscale = cscale, std.err = std.err, corr = corr, method = method,
      warn.inf = warn.inf),
    amix = fbvcamix(x = x, u = u, start = start, ..., cshape = cshape,
      cscale = cscale, std.err = std.err, corr = corr, method = method,
      warn.inf = warn.inf))
}

fbvppot <- function(x, u, model = c("log", "bilog", "alog", "neglog", "negbilog", "aneglog", "ct", "hr", "amix"), start, ..., sym = FALSE, cshape = cscale, cscale = FALSE, std.err = TRUE, corr = FALSE, method = "BFGS", warn.inf = TRUE) {
  model <- match.arg(model)
  if(model %in% c("alog","aneglog","amix"))
    stop("This model is not appropriate for poisson likelihood")
  if(sym && (model != "ct"))
    warning("Argument `sym' was ignored")
  switch(model,
    log = fbvplog(x = x, u = u, start = start, ..., cshape = cshape,
      cscale = cscale, std.err = std.err,
      corr = corr, method = method, warn.inf = warn.inf),
    bilog = fbvpbilog(x = x, u = u, start = start, ..., cshape = cshape,
      cscale = cscale, std.err = std.err, corr =
      corr, method = method, warn.inf = warn.inf),
    neglog = fbvpneglog(x = x, u = u, start = start, ..., cshape =
      cshape, cscale = cscale, std.err = std.err,
      corr = corr, method = method, warn.inf = warn.inf),
    negbilog = fbvpnegbilog(x = x, u = u, start = start, ..., cshape =
      cshape, cscale = cscale, std.err = std.err,
      corr = corr, method = method, warn.inf = warn.inf),
    ct = fbvpct(x = x, u = u, start = start, ..., sym = sym, cshape =
      cshape, cscale = cscale, std.err = std.err,
      corr = corr, method = method, warn.inf = warn.inf),
    hr = fbvphr(x = x, u = u, start = start, ..., cshape = cshape,
      cscale = cscale, std.err = std.err, corr = corr, method = method,
      warn.inf = warn.inf))
}

################## censored likelihood fitting routines #################

fbvclog <- function(x, u, start, ..., cshape = cscale, cscale = FALSE, std.err = TRUE, corr = FALSE, method = "BFGS", warn.inf = TRUE) {

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
  bvtpost.optim(x, u, opt, nm, fixed.param, std.err, corr, spx$nat, sym = FALSE, cmar = c(cscale, cshape), model = "log", warn.inf = warn.inf)
}

fbvcbilog <- function(x, u, start, ..., cshape = cscale, cscale = FALSE, std.err = TRUE, corr = FALSE, method = "BFGS", warn.inf = TRUE) {

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
  bvtpost.optim(x, u, opt, nm, fixed.param, std.err, corr, spx$nat, sym = FALSE, cmar = c(cscale, cshape), model = "bilog", warn.inf = warn.inf)
}

fbvcalog <- function(x, u, start, ..., sym = FALSE, cshape = cscale, cscale = FALSE, std.err = TRUE, corr = FALSE, method = "BFGS", warn.inf = TRUE) {
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
  bvtpost.optim(x, u, opt, nm, fixed.param, std.err, corr,
    spx$nat, sym = sym, cmar = c(cscale, cshape), model = "alog",
    warn.inf = warn.inf)
}

fbvcneglog <- function(x, u, start, ..., cshape = cscale, cscale = FALSE, std.err = TRUE, corr = FALSE, method = "BFGS", warn.inf = TRUE) {

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
  bvtpost.optim(x, u, opt, nm, fixed.param, std.err, corr, spx$nat, sym = FALSE, cmar = c(cscale, cshape), model = "neglog", warn.inf = warn.inf)
}

fbvcnegbilog <- function(x, u, start, ..., cshape = cscale, cscale = FALSE, std.err = TRUE, corr = FALSE, method = "BFGS", warn.inf = TRUE) {

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
  bvtpost.optim(x, u, opt, nm, fixed.param, std.err, corr, spx$nat, sym = FALSE, cmar = c(cscale, cshape), model = "negbilog", warn.inf = warn.inf)
}

fbvcaneglog <- function(x, u, start, ..., sym = FALSE, cshape = cscale, cscale = FALSE, std.err = TRUE, corr = FALSE, method = "BFGS", warn.inf = TRUE) {
  
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
  bvtpost.optim(x, u, opt, nm, fixed.param, std.err, corr, spx$nat, sym = sym, cmar = c(cscale, cshape), model = "aneglog", warn.inf = warn.inf)
}

fbvcct <- function(x, u, start, ..., sym = FALSE, cshape = cscale, cscale = FALSE, std.err = TRUE, corr = FALSE, method = "BFGS", warn.inf = TRUE) {

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
  bvtpost.optim(x, u, opt, nm, fixed.param, std.err, corr, spx$nat, sym = sym, cmar = c(cscale, cshape), model = "ct", warn.inf = warn.inf)
}

fbvchr <- function(x, u, start, ..., cshape = cscale, cscale = FALSE, std.err = TRUE, corr = FALSE, method = "BFGS", warn.inf = TRUE) {

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
  bvtpost.optim(x, u, opt, nm, fixed.param, std.err, corr, spx$nat, sym = FALSE, cmar = c(cscale, cshape), model = "hr", warn.inf = warn.inf)
}

fbvcamix <- function(x, u, start, ..., cshape = cscale, cscale = FALSE, std.err = TRUE, corr = FALSE, method = "BFGS", warn.inf = TRUE) {

  nllbvcamix <- function(scale1, shape1, scale2, shape2, alpha, beta) {
    if(cshape) shape2 <- shape1
    if(cscale) scale2 <- scale1
    .C("nllbvcamix", spx$x1, spx$x2, spx$nn, spx$n, spx$thdi, spx$lambda,
      alpha, beta, scale1, shape1, scale2, shape2, dns = double(1),
      PACKAGE = "evd")$dns
  }
  param <- c("scale1", "shape1")
  if(!cscale) param <- c(param, "scale2")
  if(!cshape) param <- c(param, "shape2")
  param <- c(param, "alpha", "beta")
  nmdots <- names(list(...))
  start <- bvstart.vals(x, start, NULL, NULL, nmdots, param, NULL,
    NULL, model = "amix", obj = "bvpot", u = u)
  spx <- sep.bvdata(x, obj = "bvpot", u = u) 
  nm <- names(start)
  l <- length(nm)
  fixed.param <- list(...)[nmdots %in% param]
  if(any(!(param %in% c(nm,names(fixed.param)))))
    stop("unspecified parameters")
  prind <- c(TRUE, TRUE, !cscale, !cshape, TRUE, TRUE)
  f <- formals(nllbvcamix)[prind]
  names(f) <- param
  m <- match(nm, param)
  if(any(is.na(m))) 
    stop("`start' specifies unknown arguments") 
  formals(nllbvcamix) <- c(f[m], f[-m])
  nll <- function(p, ...) nllbvcamix(p, ...)
  if(l > 1) {
    body(nll) <- parse(text = paste("nllbvcamix(", paste("p[", 1:l, "]", collapse=", "), ", ...)"))
  }
  start.arg <- c(list(p = unlist(start)), fixed.param)
  if(warn.inf && do.call("nll", start.arg) == 1e+06)
    warning("negative log-likelihood is infinite at starting values")
  opt <- optim(start, nll, hessian = std.err, ..., method = method)
  bvtpost.optim(x, u, opt, nm, fixed.param, std.err, corr, spx$nat, sym = FALSE, cmar = c(cscale, cshape), model = "amix", warn.inf = warn.inf)
}

################## Poisson likelihood fitting routines ##################

fbvplog <- function(x, u, start, ..., cshape = cscale, cscale = FALSE, std.err = TRUE, corr = FALSE, method = "BFGS", warn.inf = TRUE) {

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
  bvtpost.optim(x, u, opt, nm, fixed.param, std.err, corr,
    spx$nat, sym = FALSE, cmar = c(cscale, cshape), model = "log",
    warn.inf = warn.inf)
}

fbvpneglog <- function(x, u, start, ..., cshape = cscale, cscale = FALSE, std.err = TRUE, corr = FALSE, method = "BFGS", warn.inf = TRUE) {

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
  bvtpost.optim(x, u, opt, nm, fixed.param, std.err, corr, spx$nat, sym = FALSE, cmar = c(cscale, cshape), model = "neglog", warn.inf = warn.inf)
}

fbvpct <- function(x, u, start, ..., sym = FALSE, cshape = cscale, cscale = FALSE, std.err = TRUE, corr = FALSE, method = "BFGS", warn.inf = TRUE) {

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
  bvtpost.optim(x, u, opt, nm, fixed.param, std.err, corr, spx$nat, sym = sym, cmar = c(cscale, cshape), model = "ct", warn.inf = warn.inf)
}

fbvpbilog <- function(x, u, start, ..., cshape = cscale, cscale = FALSE, std.err = TRUE, corr = FALSE, method = "BFGS", warn.inf = TRUE) {

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
  bvtpost.optim(x, u, opt, nm, fixed.param, std.err, corr, spx$nat, sym = FALSE, cmar = c(cscale, cshape), model = "bilog", warn.inf = warn.inf)
}

fbvpnegbilog <- function(x, u, start, ..., cshape = cscale, cscale = FALSE, std.err = TRUE, corr = FALSE, method = "BFGS", warn.inf = TRUE) {

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
  bvtpost.optim(x, u, opt, nm, fixed.param, std.err, corr, spx$nat, sym = FALSE, cmar = c(cscale, cshape), model = "negbilog", warn.inf = warn.inf)
}

fbvphr <- function(x, u, start, ..., cshape = cscale, cscale = FALSE, std.err = TRUE, corr = FALSE, method = "BFGS", warn.inf = TRUE) {

  nllbvphr <- function(scale1, shape1, scale2, shape2, dep) {
    if(cshape) shape2 <- shape1
    if(cscale) scale2 <- scale1
    .C("nllbvphr", spx$x1, spx$x2, spx$nn, spx$n, spx$thdi, spx$r1, spx$r2,
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
  spx <- sep.bvdata(x, obj = "bvpot", u = u, censored = FALSE) 
  nm <- names(start)
  l <- length(nm)
  fixed.param <- list(...)[nmdots %in% param]
  if(any(!(param %in% c(nm,names(fixed.param)))))
    stop("unspecified parameters")
  prind <- c(TRUE, TRUE, !cscale, !cshape, TRUE)
  f <- formals(nllbvphr)[prind]
  names(f) <- param
  m <- match(nm, param)
  if(any(is.na(m))) 
    stop("`start' specifies unknown arguments")
  formals(nllbvphr) <- c(f[m], f[-m])
  nll <- function(p, ...) nllbvphr(p, ...)
  if(l > 1) {
    body(nll) <- parse(text = paste("nllbvphr(", paste("p[", 1:l, "]", collapse=", "), ", ...)"))
  }
  start.arg <- c(list(p = unlist(start)), fixed.param)
  if(warn.inf && do.call("nll", start.arg) == 1e+06)
    warning("negative log-likelihood is infinite at starting values")
  opt <- optim(start, nll, hessian = std.err, ..., method = method)
  bvtpost.optim(x, u, opt, nm, fixed.param, std.err, corr,
    spx$nat, sym = FALSE, cmar = c(cscale, cshape), model = "hr",
    warn.inf = warn.inf)
}

###################### post-optimisation processing #####################

bvtpost.optim <- function(x, u, opt, nm, fixed.param, std.err, corr, nat, sym, cmar, model, warn.inf) {
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
  # Dependence chi
  if(warn.inf) {
    if(model %in% c("log", "hr", "neglog")) {
      dep <- param["dep"]
      dep.sum <- 2 * (1 - abvevd(dep = dep, model = model))
    }
    if(model %in% c("alog", "aneglog")) {
      dep <- param["dep"]
      asy <- param[c("asy1", "asy2")]
      dep.sum <- 2 * (1 - abvevd(dep = dep, asy = asy, model = model))
    }
    if(model %in% c("bilog", "negbilog", "ct", "amix")) {
      alpha <- param["alpha"]
      beta <- param["beta"]
      dep.sum <- 2 * (1 - abvevd(alpha = alpha, beta = beta, model = model))
    }
  }
  else dep.sum <- NULL
  # End dependence chi
  list(estimate = opt$par, std.err = std.err, fixed = fixed, fixed2 = fixed2, param = param, deviance = 2 * opt$value, dep.summary = dep.sum, corr = corr, convergence = opt$convergence, counts = opt$counts, message = opt$message, data = x, threshold = u, n = nrow(x), nat = nat, sym = sym, cmar = cmar, model = model)
}


########################## method functions #########################


"print.bvpot" <-  function(x, digits = max(3, getOption("digits") - 3), ...) 
{
    cat("\nCall:", deparse(x$call), "\n")
    cat("Deviance:", deviance(x), "\n")
    cat("AIC:", AIC(x), "\n")
    if(!is.null(x$dep.summary)) cat("Dependence:", x$dep.summary, "\n")
    
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
    cat("\nOptimization Information\n")
    cat("  Convergence:", x$convergence, "\n")
    cat("  Function Evaluations:", x$counts["function"], "\n")
    if(!is.na(x$counts["gradient"]))
        cat("  Gradient Evaluations:", x$counts["gradient"], "\n")
    if(!is.null(x$message)) cat("  Message:", x$message, "\n")
    cat("\n")
    invisible(x)
}





