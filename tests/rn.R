##############
# FROM RNEWS #
##############

library(evd)
options(digits = 4, width = 80)
set.seed(50)

data(sealevel) ; sl <- sealevel
tt <- (1912:1992 - 1950)/100
lg <- fbvevd(sl, model = "log", nsloc1 = tt, nsloc2 = tt)
lg
lg2 <- fbvevd(sl, model = "log")
lg2
anova(lg, lg2) 
pr <- profile(lg, "dep", xmax = 1)
pr
pcint(pr)
#plot(pr)
#plot(lg)
#plot(lg, mar = 1)
#plot(lg, mar = 2)
