###################
# FROM MANUAL 1.1 #
###################

library(evd)
options(digits = 4, width = 80)
set.seed(50)

# Section: Standard Univariate Functions

rgev(6, loc = c(2,1), scale = .5, shape = 1)
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

rmvalog(3, dep = c(.6,.5,.8,.3), asy = list(.4,0,.6,c(.3,.2),c(.1,.1),c(.4,.1),c(.2,.4,.2)), d = 3)
pmvalog(c(2, 2, 2), dep = c(.6,.5,.8,.3), asy = list(.4,.0,.6,c(.3,.2),c(.1,.1),c(.4,.1),c(.2,.4,.2)), d = 3)
tmp.quant <- matrix(rep(c(1,1.5,2), 3), ncol = 3)
pmvalog(tmp.quant, dep = c(.6,.5,.8,.3), asy = list(.4,.0,.6,c(.3,.2),c(.1,.1),c(.4,.1),c(.2,.4,.2)), d = 3)
rmvalog(3, dep = c(rep(1,6),.7,.3,.8,.7,.5), asy = list(0, 0, 0, 0, c(0,0), c(0,0), c(0,0), c(0,0), c(0,0), c(0,0), c(.2,.1,.2), c(.1,.1,.2), c(.3,.4,.1), c(.2,.2,.2), c(.4,.6,.2,.5)), d = 4)
rmvalog(3, dep = c(.6,1,.8,.3), asy = list(.4,0,.6,c(.3,.2),c(0,0),c(.4,.1),c(.3,.4,.3)), d = 3)

# Section: Fitting Univariate Distributions by Maximum Likelihood

data <- rgev(1000, loc = 0.13, scale = 1.1, shape = 0.2)
fgev(data)
fgev(data, loc = 0, scale = 1)

data2 <- rext(100, qnorm, mean = 0.56, mlen = 365)
fext(data2, list(mean = 0, sd = 1), distn = "norm", mlen = 365)
fext(data2, list(scale = 1), shape = 0.5, distn = "gamma", mlen = 365)

# Section: Fitting Bivariate Distributions by Maximum Likelihood

bvdata <- rbvlog(100, dep = 0.6, mar1 = c(1.2,1.4,0), mar2 = c(1,1.6,0.1))
fbvlog(bvdata)
fbvlog(bvdata, dep = 1)
fbvlog(bvdata, shape1 = 0, shape2 = 0)
pchisq(748.4 - 728, df = 2, lower.tail = FALSE)

# Boundary Problems

fbvalog(bvdata)
fbvalog(bvdata, asy2 = 1)$estimate
fbvalog(bvdata, method = "L-BFGS-B", lower = c(rep(-Inf, 6), 0, 0, -Inf), upper = c(rep(Inf, 6), 1, 1, Inf))$estimate

# Section: Extended Example: Oxford Data

data(oxford)
oxford.fit.trend <- fgev(oxford, nsloc = (1901:1980 - 1950)/100)
oxford.fit.trend

oxford.fit <- fgev(oxford)
oxford.fit <- fgev(oxford, nsloc = (1901:1980 - 1950)/100, loctrend = 0)
oxford.fit

fgev(oxford, shape = 0)$deviance - oxford.fit$deviance
mle <- oxford.fit$estimate
as.vector(mle[1] - mle[2]/mle[3])
range(oxford)

fext(oxford, start = list(mean = 40, sd = 1), distn = "norm", mlen = 365)
fext(oxford, start = list(scale = 1, shape = 1), distn = "gamma", mlen = 365)

# Section: Extended Example: Sea Level Data

data(sealevel)
sl <- sealevel
tvec <- (1912:1992 - 1950)/100
m1 <- fbvlog(sl, nsloc1 = tvec, nsloc2 = tvec)
m2 <- fbvlog(sl, nsloc1 = tvec)
m3 <- fbvlog(sl)

pchisq(m3$deviance - m2$deviance, df = 1, lower.tail = FALSE)
pchisq(m2$deviance - m1$deviance, df = 1, lower.tail = FALSE)
pchisq(m3$deviance - m1$deviance, df = 2, lower.tail = FALSE)

mle <- m1$estimate
as.vector(mle[1] + tvec * mle[2] - mle[3]/mle[4])
as.vector(mle[5] + tvec * mle[6] - mle[7]/mle[8])

tdframe <- data.frame(linear = tvec, quad = tvec^2)
m4 <- fbvlog(sl, nsloc1 = tdframe, nsloc2 = tvec)
m4$estimate[1:3]
m4$std.err[1:3]
pchisq(m1$deviance - m4$deviance, df = 1, lower.tail = FALSE)

fbvlog(sl, nsloc1 = tvec, nsloc2 = tvec)
fbvlog(sl, nsloc1 = tdframe, nsloc2 = tvec, loc1quad = 0)

m5 <- fbvlog(sl, nsloc1 = tvec, nsloc2 = tvec, shape1 = 0, shape2 = 0)
pchisq(m5$deviance - m1$deviance, df = 2, lower.tail = FALSE)
m6 <- fbvlog(sl, nsloc1 = tvec, nsloc2 = tvec, dep = 1)
m6$deviance - m1$deviance

#fbvall(sl, nsloc1 = tvec, nsloc2 = tvec)




