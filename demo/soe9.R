########################################################
# Examples from Ch9 of "Statistics of Extremes"        #
# Read the file demos.txt for detailed information.    #
########################################################

# Graphics Settings
if(dev.cur() <= 1) get(getOption("device"))()
opar <- par(ask = interactive() && (.Device %in% c("X11", "GTK", "gnome", "windows","quartz")))

####################
# 9.1 Introduction #
####################

# Figure 9.1 (a-b)
data(lossalae); nn <- nrow(lossalae)
lossalae <- lossalae/100000; lts <- c(1e-04, 100)
plot(lossalae, log = "xy", xlim = lts, ylim = lts)
ula <- apply(lossalae, 2, rank)/(nn + 1)
plot(ula)

#########################
# 9.2 Parametric Models #
#########################

# Figure 9.2 (a-b)
abvevd(dep = 0.5, asy = c(1,1), model = "alog", plot = TRUE)
abvevd(dep = 0.5, asy = c(0.6,0.9), model = "alog", rev = TRUE, add = TRUE, lty = 2)
abvevd(dep = 0.5, asy = c(0.8,0.5), model = "alog", rev = TRUE, add = TRUE, lty = 3)
abvevd(dep = -1/(-2), model = "neglog", plot = TRUE)
abvevd(dep = -1/(-1), model = "neglog", add = TRUE, lty = 2)
abvevd(dep = -1/(-0.5), model = "neglog", add = TRUE, lty = 3)

# Figure 9.3 (a-b)
abvevd(alpha = 1, beta = -0.2, model = "amix", plot = TRUE)
abvevd(alpha = 0.6, beta = 0.1, model = "amix", add = TRUE, lty = 2)
abvevd(alpha = 0.2, beta = 0.2, model = "amix", add = TRUE, lty = 3)
abvevd(dep = 1/1.25, model = "hr", plot = TRUE)
abvevd(dep = 1/0.83, model = "hr", add = TRUE, lty = 2)
abvevd(dep = 1/0.5, model = "hr", add = TRUE, lty = 3)

#############################
# 9.3 Component-wise Maxima #
#############################

# Random Cwise Maxima Derivation 
set.seed(50)
cmla <- lossalae[sample(nn),]
xx <- rep(1:50, each = 30)
cmla <- cbind(tapply(cmla[,1], xx, max), tapply(cmla[,2], xx, max))
colnames(cmla) <- colnames(lossalae)

# Figure 9.4 (a-b)
plot(lossalae, log = "xy", xlim = lts, ylim = lts, col = "grey")
points(cmla)
ecmla <- -log(apply(cmla,2,rank)/51)
plot(ecmla)

# Figure 9.5 (a-b)
# Using Empirical Margins For Figure 9.5(a)
# Using Gev Margins For Figure 9.5(b)
abvnonpar(data = cmla, epmar = TRUE, method = "pick", rev = TRUE, plot = TRUE, lty = 3)
abvnonpar(data = cmla, epmar = TRUE, method = "pick", rev = TRUE, add = TRUE, madj = 1, lty = 2)
abvnonpar(data = cmla, epmar = TRUE, method = "pick", rev = TRUE, add = TRUE, madj = 2, lty = 4)
abvnonpar(data = cmla, epmar = TRUE, rev = TRUE, add = TRUE, lty = 1)
m1 <- fbvevd(cmla, asy1 = 1, model = "alog")
m2 <- fbvevd(cmla, model = "log")
m3 <- fbvevd(cmla, model = "bilog")
plot(m1, which = 4, rev = TRUE, nplty = 3)
plot(m2, which = 4, rev = TRUE, nplty = 3, lty = 2, add = TRUE)
plot(m3, which = 4, rev = TRUE, nplty = 3, lty = 4, add = TRUE)

# Table 9.1 (Using Gev Margins)
fitted(m2); std.errors(m2); deviance(m2)/2
fitted(m1); std.errors(m1); deviance(m1)/2
fitted(m3); std.errors(m3); deviance(m3)/2

# Tawn Score Test (Using Empirical Margins)
rsm <- rowSums(ecmla); rsl <- rowSums(ecmla * log(ecmla))
tawn <- rsl - log(apply(ecmla, 1, prod)) - (rsm - 2) * log(rsm) - 1/rsm
pnorm((25 * log(50))^(-1/2) * sum(tawn))

# Tawn Score Test (Using Gev Margins)
mmles <- list(fitted(fgev(cmla[,1])), fitted(fgev(cmla[,2])))
ecmla <- mtransform(cmla, mmles)
rsm <- rowSums(ecmla); rsl <- rowSums(ecmla * log(ecmla))
tawn <- rsl - log(apply(ecmla, 1, prod)) - (rsm - 2) * log(rsm) - 1/rsm
pnorm((25 * log(50))^(-1/2) * sum(tawn))

# Likelihood Ratio Tests (Using Gev Margins)
anova(m3, m2)
anova(m1, m2, half = TRUE)

# Figure 9.6 (Using Empirical Margins)
lts <- c(0.01,100)
plot(lossalae, log = "xy", col = "grey", xlim = lts, ylim = lts)
points(cmla)
qcbvnonpar(c(0.98,0.99,0.995), data = cmla, epmar = TRUE, mint = 30, add = TRUE)

# Figure 9.6 (Using Gev Margins)
plot(lossalae, log = "xy", col = "grey", xlim = lts, ylim = lts)
points(cmla)
qcbvnonpar(c(0.98,0.99,0.995), data = cmla, mint = 30, add = TRUE)

#################################
# 9.4 Excesses Over a Threshold #
#################################

# Figure 9.7 (a-b)
fla <- -1/log(ula); pla <- 1/(1-ula)
rr <- rowSums(fla); ww <- fla/rr
rro <- sort(rr, decreasing = TRUE)[-1]
k <- 1:(nn-1)
plot(k, rro*k/nn, ylab = "")
abline(h = 2, v = 337)
xx <- yy <- seq(0, 1, len = 100)
for(k in 1:100) yy[k] <- sum(rr > rro[337] & ww[,1] <= xx[k])
plot(xx, 2/337 * yy, type = "l", ylab = "H([0,w])")
abline(h = c(0,2))

# Marginal Fits
thresh <- apply(lossalae, 2, sort, decreasing = TRUE)[170,]
fitted(fpot(lossalae[,1], thresh[1]))
fitted(fpot(lossalae[,2], thresh[2]))

# Table 9.2
m1 <- fbvpot(lossalae, thresh, model = "alog", asy1 = 1)
m2 <- fbvpot(lossalae, thresh, model = "bilog")
m3 <- fbvpot(lossalae, thresh, model = "bilog", likelihood = "poisson")
fitted(m1); std.errors(m1)
fitted(m2); std.errors(m2)
fitted(m3); std.errors(m3)

# Figure 9.8 (a-b)
abvnonpar(data = lossalae, method = "pot", k = 337, epmar = TRUE, plot = TRUE, rev = TRUE, lty = 3)
plot(m1, which = 2, rev = TRUE, add = TRUE)
plot(m2, which = 2, rev = TRUE, add = TRUE, lty = 4)
plot(m3, which = 2, rev = TRUE, add = TRUE, lty = 2)
lts <- c(1e-04, 100)
plot(lossalae, log = "xy", col = "grey", xlim = lts, ylim = lts)
plot(m1, which = 3, p = c(0.98,0.99,0.995,0.999), tlty = 0, add = TRUE)

# Figure 9.9
chiplot(lossalae, ylim1 = c(-1,1), nq = 200, qlim = c(0.02,0.98), which = 1, xlab = "u", ylab1 = "Chi(u)", main1 = "", spcases = TRUE)
chiplot(lossalae, nq = 200, qlim = c(0.02,0.98), which = 2, xlab = "u", ylab2 = "Chib(u)", main2 = "", spcases = TRUE)

# Figure 9.10
fla <- apply(fla, 1, min); pla <- apply(pla, 1, min)
thresh <- quantile(fla, probs = c(0.025, 0.975))
tcplot(fla, thresh, nt = 100, pscale = TRUE, which = 2, vci = FALSE, cilty = 2, type = "l", ylim = c(-0.2,1.2), ylab = "eta")
abline(h = c(0,1))
thresh <- quantile(pla, probs = c(0.025, 0.975))
tcplot(pla, thresh, nt = 100, pscale = TRUE, which = 2, vci = FALSE, cilty = 2, type = "l", ylim = c(-0.2,1.2), ylab = "eta")
abline(h = c(0,1))

# Example of Likelihood Ratio Test
thresh <- quantile(fla, probs = 0.8)
m1 <- fpot(fla, thresh = thresh)
m2 <- fpot(fla, thresh = thresh, shape = 1)
anova(m1, m2, half = TRUE)

# Figure 9.11 Omitted
# Return Graphics Settings
par(opar)
