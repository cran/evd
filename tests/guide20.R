###################
# FROM MANUAL 2.0 #
###################

library(evd)
options(digits = 4, width = 80)
set.seed(50)

# Section: Univariate Distributions

rgev(6, loc = c(20,1), scale = .5, shape = 1)
qrweibull(seq(0.1, 0.4, 0.1), 2, 0.5, 1, lower.tail = FALSE)
qrweibull(seq(0.9, 0.6, -0.1), loc = 2, scale = 0.5, shape = 1)
pfrechet(2:6, 2, 0.5, 1)
pfrechet(2:6, 2, 0.5, 1, low = FALSE)
drweibull(-1:3, 2, 0.5, log = TRUE)
dgumbel(-1:3, 0, 1)

rextreme(1, distn = "norm", sd = 2, mlen = 20, largest = FALSE)
min(rnorm(20, mean = 0, sd = 2))
rextreme(4, distn = "exp", rate = 1, mlen = 5)
rextreme(4, distn = "exp", mlen = 5)
pextreme(c(.4, .5), distn = "norm", mean = 0.5, sd = c(1, 2), mlen = 4)
dextreme(c(1, 4), distn = "gamma", shape = 1, scale = 0.3, mlen = 100)

rorder(1, distn = "norm", mlen = 20, j = 2)
porder(c(1, 2), distn = "gamma", shape = c(.5, .7), mlen = 10, j = 2)
dorder(c(1, 2), distn = "gamma", shape = c(.5, .7), mlen = 10, j = 2)

# Section: Bivariate Extreme Value Distributions

rbvevd(3, dep = .8, asy = c(.4, 1), model = "alog")
rbvevd(3, alpha = .5, beta = 1.2, model = "negb", mar1 = rep(1, 3))
pbvevd(c(1, 1.2), dep = .4, asy = c(.4, .6), model = "an", mar1 = rep(1, 3))
tmp.quant <- matrix(c(1,1.2,1,2),ncol = 2, byrow = TRUE)
tmp.mar <- matrix(c(1,1,1,1.2,1.2,1.2), ncol = 3, byrow = TRUE)
pbvevd(tmp.quant, dep = .4, asy = c(.4, .6), model = "an", mar1 = tmp.mar)
dbvevd(c(1, 1.2), alpha = .2, beta = .6, model = "ct", mar1 = rep(1, 3))
dbvevd(tmp.quant, alpha = 0.2, beta = 0.6, model = "ct", mar1 = tmp.mar)

# Section: Multivariate Extreme Value Distributions

rmvevd(3, dep = .6, model = "log", d = 5)
tmp.mar <- matrix(c(1,1,1,1,1,1.5,1,1,2), ncol = 3, byrow = TRUE)
rmvevd(3, dep = .6, d = 5, mar = tmp.mar)
tmp.quant <- matrix(rep(c(1,1.5,2), 5), ncol = 5)
pmvevd(tmp.quant, dep = .6, d = 5, mar = tmp.mar)
dmvevd(tmp.quant, dep = .6, d = 5, mar = tmp.mar, log = TRUE)

asy <- list(.4, 0, .6, c(.3,.2), c(.1,.1), c(.4,.1), c(.2,.4,.2))
rmvevd(3, dep = c(.6,.5,.8,.3), asy = asy, model = "alog", d = 3)
dmvevd(c(2, 2, 2), dep = c(.6,.5,.8,.3), asy = asy, model = "a", d = 3)
tmp.quant <- matrix(rep(c(1,1.5,2), 3), ncol = 3)
pmvevd(tmp.quant, dep = c(.6,.5,.8,.3), asy = asy, model = "a", d = 3)
asy <- list(0, 0, 0, 0, c(0,0), c(0,0), c(0,0), c(0,0), c(0,0), c(0,0), c(.2,.1,.2), c(.1,.1,.2), c(.3,.4,.1), c(.2,.2,.2), c(.4,.6,.2,.5))
rmvevd(3, dep = c(rep(1,6),.7,.3,.8,.7,.5), asy = asy, model = "alog", d = 4)

# Section: Dependence Functions

bvlsm <- rmvevd(100, dep = 0.6, model = "log", d = 2)
tvlsm <- rmvevd(100, dep = 0.6, model = "log", d = 3)
abvpar(seq(0,1,0.25), dep = 0.3, asy = c(.7,.9), model = "alog")
abvnonpar(seq(0,1,0.25), data = bvlsm)

# Section: Stochastic Processes

marma(100, p = 1, q = 1, psi = 0.75, theta = 0.65)
mar(100, psi = 0.85, n.start = 20)
mma(100, q = 2, theta = c(0.75, 0.8))
evmc(100, alpha = 0.1, beta = 0.1, model = "bilog")
evmc(100, dep = 10, model = "hr", margins = "exp")

# Section: Fitting Univariate Distributions

data1 <- rgev(1000, loc = 0.13, scale = 1.1, shape = 0.2)
m1 <- fgev(data1)
m1
m2 <- fgev(data1, loc = 0, scale = 1)
m2
fitted(m1)
std.errors(m1)
anova(m1,m2)

d2 <- rextreme(100, distn = "norm", mean = 0.56, mlen = 365)
sv <- list(mean = 0, sd = 1)
nm <- fextreme(d2, start = sv, distn = "norm", mlen = 365)
fitted(nm)
d3 <- rorder(100, distn = "norm", mean = 0.56, mlen = 365, j = 2)
sv <- list(mean = 0, sd = 1)
nm2 <- forder(d3, sv, distn = "norm", mlen = 365, j = 2)
fitted(nm2)

# Section: Fitting Bivariate Extreme Value Distributions

bvdata <- rbvlog(100, dep = 0.6, mar1 = c(1.2,1.4,0), mar2 = c(1,1.6,0.1))
m1 <- fbvevd(bvdata, model = "log")
m1
m2 <- fbvevd(bvdata, model = "log", dep = 1)
fitted(m2)
std.errors(m2)
deviance(m2)
m3 <- fbvevd(bvdata, model = "log", shape1 = 0, shape2 = 0)
anova(m1, m3)
m4 <- fbvevd(bvdata, model = "alog")
fitted(m4)
mb <- fbvevd(bvdata, model = "alog", asy2 = 1)
round(fitted(mb), 3)
up <- c(rep(Inf, 6), 1, 1, 1)
mb <- fbvevd(bvdata, model = "alog", method = "L-BFGS-B", upper = up)
round(fitted(mb), 3)

# Section: Example: Oxford Temperature Data

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

ox.qfit <- fgev(ox, prob = 0.1)
ox.qprof <- profile(ox.qfit, which = "quantile")
ox.qfit <- fgev(ox, prob = 0)
#ox.qprof <- profile(ox.qfit, which = "quantile", conf = 0.99)

# Section: Example: Sea Level Data

data(sealevel) ; sl <- sealevel
tt <- (1912:1992 - 1950)/100
m1 <- fbvevd(sl, model = "log", nsloc1 = tt, nsloc2 = tt)
m2 <- fbvevd(sl, model = "log", nsloc1 = tt)
m3 <- fbvevd(sl, model = "log")
anova(m1, m2, m3)

tdframe <- data.frame(trend = tt, quad = tt^2)
m4 <- fbvevd(sl, model = "log", nsloc1 = tdframe, nsloc2 = tt)
m5 <- fbvevd(sl, model = "log", nsloc1 = tt, nsloc2 = tdframe)
m6 <- fbvevd(sl, model = "log", nsloc1 = tdframe, nsloc2 = tdframe)

m7 <- fbvevd(sl, model = "log", nsloc1 = tt, nsloc2 = tt, dep = 1)
anova(m1, m7)
m8 <- fbvevd(sl, model = "log", nsloc1 = tt, nsloc2 = tt, shape1 = 0, shape2 = 0)
anova(m1, m8)
#m1.prof <- profile(m1, which = "dep", xmax = 1)



