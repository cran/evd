###################
# FROM MANUAL 1.0 #
###################

library(evd)
options(digits = 4, width = 80)
set.seed(50)

# Section: Standard Univariate Functions

rrweibull(6, loc = 2, scale = .5, shape = c(1, 1.2))
qrweibull(seq(0.1, 0.4, 0.1), 2, 0.5, 1, lower.tail = FALSE)
qrweibull(seq(0.9, 0.6, -0.1), loc = 2, scale = 0.5, shape = 1)
prweibull(-1:3, 2, 0.5, 1)
prweibull(-1:3, 2, 0.5, 1, low = FALSE)
drweibull(-1:3, loc = 2, scale = 0.5, shape = 1)
drweibull(-1:3, 2, 0.5, 1, log = TRUE)

rext(4, qexp, rate = 1, mlen = 5)
rext(4, distn = "exp", rate = 1, mlen = 5)
rext(4, distn = "exp", mlen = 5)
rext(1, distn = "norm", mean = 0.5, sd = 2, mlen = 20)
max(rnorm(20, 0.5, 2))
rext(1, distn="norm", sd = 2, mlen = 20, largest = FALSE)
min(rnorm(20, 0, 2))
pext(c(.4, .5), distn="norm", sd = c(1, 2), mlen = 4)
pext(c(.4, .5), distn="norm", mean = 0, sd = c(1, 2), mlen = 4)
dext(c(1, 4), distn="gamma", shape = 1, scale = 0.3, mlen = 100)

rorder(1, distn = "norm", mlen = 20, j = 2)
rorder(1, distn = "norm", mlen = 20, j = 19, largest = FALSE)
porder(c(1, 2), distn="gamma", shape =c(.5, .7), mlen = 10, j = 2)
dorder(c(1, 2), distn="gamma", shape =c(.5, .7), mlen = 10, j = 2)

# Section: Standard Bivariate and Multivariate Functions

rbvalog(3, dep = .8, asy = c(.4, 1))
rbvalog(3, dep = .8, asy = c(.4, 1), mar1 = c(1, 1, 1))
pbvalog(c(1, 1.2), dep = .4, asy = c(.4, .6), mar1 = c(1, 1, 1))
tmp.quant <- matrix(c(1,1.2,1,2),ncol = 2, byrow = TRUE)
tmp.mar <- matrix(c(1,1,1,1.2,1.2,1.2), ncol = 3, byrow = TRUE)
pbvalog(tmp.quant, dep = .4, asy = c(.4, .6), mar1 = tmp.mar)
dbvalog(c(1, 1.2), dep = .4, asy = c(.4, .6), mar1 = c(1, 1, 1))
dbvalog(tmp.quant, dep = .4, asy = c(.4, .6), mar1 = tmp.mar)
abvalog(dep = .3, asy= c(.7, .9))
abvalog(seq(0, 1, 0.25), dep = .3, asy = c(.7, .9))

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

rmvalog(3, dep = c(.6,1,.8,.3), asy = list(.4,0,.6,c(.3,.2),c(.0,.0),c(.4,.1),c(.3,.4,.3)), d = 3)

# Section: Fitting Distributions by Maximum Likelihood

# Subsection: Univariate Fitting

data <- rgev(1000, loc = 0.13, scale = 1.1, shape = 0.2)
fgev(data, start = list(loc = 0, scale = 1, shape = 0), method = "BFGS")
fgev(data, start = list(loc = 0, scale = 1), shape = 0, method = "BFGS")
fgev(data, start = list(loc = 0), scale = 1, shape = 0)
fgev(data, start = list(loc = 0))

-2*sum(dgev(data, loc = 0.1594, scale = 1.1422, shape = 0.2101, log = TRUE))

somedata <- rgev(100,1,1,0.05)
fgumbel(somedata, start = list(loc = 0.5, scale = 2))
fgev(somedata, start = list(loc = 0.5, scale = 2))

ffrechet(data, start = list(loc=-2,scale=1,shape=5), method ="BFGS")
data2 <- rgev(1000, loc = 0.13, scale = 1.1, shape = -0.25)
# CHANGE: starting values
frweibull(data2, start = list(loc = 5, scale = 5, shape = 10), method="BFGS")
fgev(data2, start = list(loc=0.13,scale=1.1,shape=0), method="BFGS")

data3 <- rext(100, qnorm, mean = 0.56, mlen = 365)
fext(data3, list(mean = 0, sd = 1), distn = "norm", mlen = 365)
fext(data3, list(rate = 1), distn = "exp", mlen = 365)
fext(data3, list(scale = 1), shape = 0.5, distn = "gamma", mlen = 365)

# Subsection: Bivariate Fitting

bvdata <- rbvlog(100, dep = 0.6, mar1 = c(1.2,1.4,0.4), mar2 = c(1.2,1.4,0.4))
fbvlog(bvdata, start = list(mar1 = c(2,1,0), mar2 = c(1,1,0), dep = 0.75), control = list(maxit = 2000))
fbvlog(bvdata, start = list(loc1=2, scale1=1, shape1=0, loc2=1, scale2=1, shape2=0, dep=0.75), control = list(maxit = 2000))
fbvlog(bvdata, start = list(mar1 = c(2,1,0), mar2 = c(1,1,0), dep = 0.75), method = "BFGS")$counts

fbvlog(bvdata, start = list(mar1 = c(2,1,0), mar2 = c(1,1,0)), dep = 1, method="BFGS")
fgev(bvdata[,1], start = list(loc=1,scale=1,shape=0), method="BFGS")$estimate
fgev(bvdata[,2], start = list(loc=1,scale=1,shape=0), method="BFGS")$estimate
 
fbvlog(bvdata, start = list(loc1=1, scale1=1, loc2=1, scale2=1, dep=0.75), shape1 = 0, shape2 = 0, method = "BFGS")
fbvlog(bvdata, start = list(mar2 = c(1,1,0), dep = 0.75), loc1 = 1, scale1 = 1.5, shape1 = 0.5, method = "BFGS")

fbvalog(bvdata, start = list(mar1 = c(1,1,0), mar2 = c(1,1,0), asy = c(.7,.7), dep = 0.75), method = "BFGS")
fbvalog(bvdata, start = list(mar1 = c(0.9,1.2,0.2), mar2 = c(1,1.2,0.3), asy2 = .72, dep = 0.58), asy1 = 1, method = "BFGS")
# CHANGE: std.err = FALSE
fbvalog(bvdata, start = list(mar1 = c(0.9,1.2,0.2), mar2 = c(1,1.2,0.3), asy = c(.99,.72), dep = 0.58), method = "L-BFGS-B", lower = c(rep(-Inf, 6), 0, 0, -Inf), upper = c(rep(Inf, 6), 1, 1, 1), std.err = FALSE)
fbvalog(bvdata, start = list(mar1 = c(1.5,1.4,0.1), mar2 = c(1.4,1.4,0.2), 
dep = 0.73), asy1 = 1, asy2 = 1, method="BFGS")

# Subsection: A Univariate Example

data(oxford)
sqrt(6 * var(oxford))/pi
mean(oxford) - 0.58 * sqrt(6 * var(oxford))/pi
oxford.fit <- fgev(oxford, start = list(loc=83.5, scale=3.5, shape=0))
oxford.fit
fgev(oxford, start = list(loc=83.8, scale=4.25))$deviance - oxford.fit$deviance

mle <- oxford.fit$estimate
as.vector(mle[1] - mle[2]/mle[3])
range(oxford)

fext(oxford, start = list(mean = 40, sd = 1), distn = "norm", mlen = 365)
fext(oxford, start = list(scale = 1, shape = 1), distn = "gamma", mlen = 365)

# Subsection: A Bivariate Example

data(sealevel)
sl <- sealevel
# CHANGE: na.rm = TRUE
sqrt(6 * c(var(sl[,1], na.rm = TRUE), var(sl[,2], na.rm = TRUE)))/pi
c(mean(sl[,1], na.rm = TRUE), mean(sl[,2], na.rm = TRUE)) - 0.58 * c(0.21, 0.24)

tmp <- fbvlog(sl, start = list(mar1 = c(3.6, 0.2, 0), mar2 = c(2.6, 0.25, 0)), dep = 1, method = "BFGS")
tmp <- fbvalog(sl, start = list(mar1 = c(3.6, 0.2, 0), mar2 = c(2.6, 0.25, 0)), asy1 = 1, asy2 = 1, dep = 1, method = "BFGS")
# CHANGE: dependence parameter fixed at lower limit
tmp <- fbvhr(sl, start = list(mar1 = c(3.6, 0.2, 0), mar2 = c(2.6, 0.25, 0)), dep = 0.2, method = "BFGS")
tmp$estimate

sl.fit <- fbvalog(sl, start = list(mar1 = c(3.6, 0.2, 0), mar2 = c(2.6, 0.25, 0), asy = c(0.8, 0.8), dep = 0.6), method = "BFGS", control = list(trace=1))
sl.fit

fbvalog(sl, start = list(loc1 = 3.6, scale1 = 0.19, loc2 = 2.6, scale2 = 0.2, asy = c(0.7, 0.45), dep = 0.24), method = "BFGS")$deviance - sl.fit$deviance
tmp$deviance - sl.fit$deviance

# CHANGE: na.rm = TRUE
mle <- sl.fit$estimate
as.vector(mle[1] - mle[2]/mle[3])
range(sl[,1], na.rm = TRUE)
as.vector(mle[4] - mle[5]/mle[6])
range(sl[,2], na.rm = TRUE)

sl.fit2 <- fbvalog(sl, start = list(mar1 = c(3.6, 0.19, -0.04), mar2 = c(2.6, 0.2, 0.09), dep = 0.6), asy1 = 1, asy2 = 1, method = "BFGS")
sl.fit2 <- fbvlog(sl, start = list(mar1 = c(3.6, 0.19, -0.04), mar2 = c(2.6, 0.2, 0.09), dep = 0.6), method = "BFGS")
sl.fit2$estimate
sl.fit2$deviance
sl.fit2$deviance - sl.fit$deviance
abvalog(dep = 0.24316, asy = c(0.69554, 0.44967))
abvlog(dep = 0.62474)

fbvhr(sl, start = list(mar1 = c(3.6, 0.19, -0.04), mar2 = c(2.6, 0.2, 0.09),  dep = 1), method = "BFGS", control = list(trace=1))
abvhr(dep = 1.253410)
fbvaneglog(sl, start = list(mar1 = c(3.6, 0.19, -0.04), mar2 = c(2.6, 0.2, 0.09), dep = 1, asy = c(0.8,0.8)), method = "BFGS", control = list(trace=1))
abvaneglog(dep = 3.44762, asy = c(0.69796,0.44601))
fbvaneglog(sl, start = list(mar1 = c(3.6, 0.19, -0.04), mar2 = c(2.6, 0.2, 0.09), dep = 1), asy1 = 1, asy2 = 1, method = "BFGS", control = list(trace=1))
abvneglog(dep = 0.87394)

fbvalog(sl, start = list(mar1 = c(3.63, 0.18, -0.04), mar2 = c(2.63, 0.2, 0.09), dep = 0.24), asy1 = 0.4492, asy2 = 0.4492, method = "BFGS")









