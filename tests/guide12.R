###################
# FROM MANUAL 1.2 #
###################

library(evd)
options(digits = 4, width = 80)
set.seed(50)

# Section: Standard Univariate Functions

rgev(6, loc = c(20,1), scale = .5, shape = 1)
qrweibull(seq(0.1, 0.4, 0.1), 2, 0.5, 1, lower.tail = FALSE)
qrweibull(seq(0.9, 0.6, -0.1), loc = 2, scale = 0.5, shape = 1)
pfrechet(2:6, 2, 0.5, 1)
pfrechet(2:6, 2, 0.5, 1, low = FALSE)
drweibull(-1:3, 2, 0.5, log = TRUE)
dgumbel(-1:3, 0, 1)

rext(4, qexp, rate = 1, mlen = 5)
rext(4, distn = "exp", rate = 1, mlen = 5)
rext(4, distn = "exp", mlen = 5)
rext(1, distn = "norm", sd = 2, mlen = 20, largest = FALSE)
min(rnorm(20, mean = 0, sd = 2))
pext(c(.4, .5), distn = "norm", mean = 0.5, sd = c(1, 2), mlen = 4)
dext(c(1, 4), distn = "gamma", shape = 1, scale = 0.3, mlen = 100)

rorder(1, distn = "norm", mlen = 20, j = 2)
rorder(1, distn = "norm", mlen = 20, j = 19, largest = FALSE)
porder(c(1, 2), distn = "gamma", shape = c(.5, .7), mlen = 10, j = 2)
dorder(c(1, 2), distn = "gamma", shape = c(.5, .7), mlen = 10, j = 2)

# Section: Standard Bivariate Functions

rbvalog(3, dep = .8, asy = c(.4, 1))
rbvnegbilog(3, alpha = .5, beta = 1.2, mar1 = c(1, 1, 1))
pbvaneglog(c(1, 1.2), dep = .4, asy = c(.4, .6), mar1 = c(1, 1, 1))
tmp.quant <- matrix(c(1,1.2,1,2),ncol = 2, byrow = TRUE)
tmp.mar <- matrix(c(1,1,1,1.2,1.2,1.2), ncol = 3, byrow = TRUE)
pbvaneglog(tmp.quant, dep = .4, asy = c(.4, .6), mar1 = tmp.mar)
dbvct(c(1, 1.2), alpha = .2, beta = .6, mar1 = c(1, 1, 1))
dbvct(tmp.quant, alpha = 0.2, beta = 0.6, mar1 = tmp.mar)
abvlog(dep = .3)
abvlog(seq(0, 1, 0.25), dep = .3)

# Section: Standard Multivariate Functions

rmvlog(3, dep = .6, d = 5)
tmp.mar <- matrix(c(1,1,1,1,1,1.5,1,1,2), ncol = 3, byrow = TRUE)
rmvlog(3, dep = .6, d = 5, mar = tmp.mar)
tmp.quant <- matrix(rep(c(1,1.5,2), 5), ncol = 5)
pmvlog(tmp.quant, dep = .6, d = 5, mar = tmp.mar)

asy <- list(.4, 0, .6, c(.3,.2), c(.1,.1), c(.4,.1), c(.2,.4,.2))
rmvalog(3, dep = c(.6,.5,.8,.3), asy = asy, d = 3)
pmvalog(c(2, 2, 2), dep = c(.6,.5,.8,.3), asy = asy, d = 3)
tmp.quant <- matrix(rep(c(1,1.5,2), 3), ncol = 3)
pmvalog(tmp.quant, dep = c(.6,.5,.8,.3), asy = asy, d = 3)
asy <- list(0, 0, 0, 0, c(0,0), c(0,0), c(0,0), c(0,0), c(0,0), c(0,0), c(.2,.1,.2), c(.1,.1,.2), c(.3,.4,.1), c(.2,.2,.2), c(.4,.6,.2,.5))
rmvalog(3, dep = c(rep(1,6),.7,.3,.8,.7,.5), asy = asy, d = 4)
asy <- list(.4, 0, .6, c(.3,.2), c(0,0), c(.4,.1), c(.3,.4,.3))
rmvalog(3, dep = c(.6,1,.8,.3), asy = asy, d = 3)

# Section: Fitting Univariate Distributions by Maximum Likelihood

data1 <- rgev(1000, loc = 0.13, scale = 1.1, shape = 0.2)
m1 <- fgev(data1)
m1
m2 <- fgev(data1, loc = 0, scale = 1)
m2
fitted(m1)
std.errors(m1)
anova(m1,m2)

data2 <- rext(100, qnorm, mean = 0.56, mlen = 365)
nm <- fext(data2, list(mean = 0, sd = 1), distn = "norm", mlen = 365)
fitted(nm)
ga <- fext(data2, list(scale = 1), shape = 0.5, distn = "gamma", mlen = 365, method="L-BFGS-B", lower = 0.01)
fitted(ga)

# Section: Fitting Bivariate Distributions by Maximum Likelihood

bvdata <- rbvlog(100, dep = 0.6, mar1 = c(1.2,1.4,0), mar2 = c(1,1.6,0.1))
m1 <- fbvlog(bvdata)
m1
m2 <- fbvlog(bvdata, dep = 1)
fitted(m2)
m3 <- fbvlog(bvdata, shape1 = 0, shape2 = 0)
anova(m1, m3)

# Boundary Problems

m4 <- fbvalog(bvdata)
fitted(m4)
fitted(fbvalog(bvdata, asy1 = 1))
upper <- c(rep(Inf, 6), 1, 1, 1) 
fitted(fbvalog(bvdata, method = "L-BFGS-B", upper = upper))

# Section: Extended Example: Oxford Data

data(oxford) ; ox <- oxford
ox.fit <- fgev(ox)
tt <- (1901:1980 - 1950)/100
ox.fit.trend <- fgev(ox, nsloc = tt)
fitted(ox.fit.trend)
std.errors(ox.fit.trend)
ox.fit

ox.fit.gum <- fgev(ox, shape = 0)
anova(ox.fit, ox.fit.gum)
ox.prof <- profile(ox.fit)
#ox.prof2d <- profile2d(ox.fit, ox.prof, which = c("scale", "shape"))

ox.qfit <- fgev.quantile(ox, prob = 0.1)
ox.qprof <- profile(ox.qfit, which = "quantile")

ox.nm <- fext(ox, list(mean = 40, sd = 1), distn = "norm", mlen = 365)
fitted(ox.nm)
ox.ga <- fext(ox, list(scale = 1, shape = 1), distn = "gamma", mlen = 365)
fitted(ox.ga)

# Section: Extended Example: Sea Level Data

data(sealevel) ; sl <- sealevel
tt <- (1912:1992 - 1950)/100
m1 <- fbvlog(sl, nsloc1 = tt, nsloc2 = tt)
m2 <- fbvlog(sl, nsloc1 = tt)
m3 <- fbvlog(sl)
anova(m1, m2, m3)

tdframe <- data.frame(trend = tt, quad = tt^2)
m4 <- fbvlog(sl, nsloc1 = tdframe, nsloc2 = tt)
m5 <- fbvlog(sl, nsloc1 = tt, nsloc2 = tdframe)
m6 <- fbvlog(sl, nsloc1 = tdframe, nsloc2 = tdframe)

m7 <- fbvlog(sl, nsloc1 = tt, nsloc2 = tt, dep = 1)
deviance(m7) - deviance(m1)
m8 <- fbvlog(sl, nsloc1 = tt, nsloc2 = tt, shape1 = 0, shape2 = 0)
anova(m1, m8)

#m1.prof <- profile(m1, which = "dep", xmax = 1)
#m.all <- fbvall(sl, nsloc1 = tt, nsloc2 = tt)
#fitted(m.all)
#m.all$dep.summary
#m.all$criteria

m9 <- fbvalog(sl, nsloc1 = tt, nsloc2 = tt)
