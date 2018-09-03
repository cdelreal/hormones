# calculate standard curve

ss <- read.table("raw.csv", header = TRUE, sep = ',')

# removing extraneous standards from the curve
# killna <- which(ss$standard_conc %in% c(50))
# ss$standard_conc[killna] <- NA

# is blank already subtracted?
deseblank <- which(ss$lab == "blank")
blank <- mean(ss[deseblank,]$rawod)
# od <- ss$rawod - blank
od <- ss$rawod

desensb <- which(ss$lab == 'nsb')
nsb <- mean(od[desensb])
deseb0 <- which(ss$lab == 'b0')
#subtract nsb from b0

b0 <- mean(od[deseb0])
netod <- od

if(length(desensb) > 0) {
	b0 <- b0 - nsb
	netod <- od - nsb #do i need to subtract nonspecific binding?
}

pbb0 <- netod / b0

desestands <- which(!is.na(ss$standard_conc))
standards <- netod[desestands]
standardlabels <- as.character(ss$lab[desestands])
conc <- ss$standard_conc[desestands]
conc <- as.vector(tapply(conc, standardlabels, mean, na.rm = TRUE))
standards_mean <- as.vector(tapply(standards, standardlabels, mean, na.rm = TRUE))

# plot(conc, standards_mean, log = 'x')

library(nplr)

# do it by percent binding
standards_pbb0 <- standards_mean / b0

mod_percentbinding <- nplr(conc, standards_pbb0, npars = 4)
ypred <- getEstimates(mod_percentbinding, standards_pbb0)$x
yexp <- conc
plot(yexp, ypred, log = 'xy')
abline(0,1)
cbind(ypred, yexp, round(100*ypred / yexp))


# calculate concentration from OD
# desesamples <- which(substring(ss$lab, 1, 1) == "b" | substring(ss$lab, 1, 1) == "m")

desesamples <- which(ss$dilution == 1)

absconc_pgmL <- getEstimates(mod_percentbinding, pbb0[desesamples])$x
# par(mfrow = c(1,2))
plot(c(conc, 10^mod_percentbinding@xCurve), c(standards_pbb0, mod_percentbinding@yCurve), type = 'n', log = 'x', las = 1)
lines(10^mod_percentbinding@xCurve, mod_percentbinding@yCurve)
points(conc, standards_pbb0, pch = 16)
# gof <- attr(mod_percentbinding, "goodness")
# title(paste("gof = ", round(gof, 3), sep = ""))


# par(oma = c(0, 0, 0, 2))
plot(c(conc, absconc_pgmL), c(standards_pbb0, pbb0[desesamples]), log = 'x', type = 'n', las = 1, xlab = "sample concentration (pg/mL)", ylab = "percent binding")
points(conc, standards_pbb0, pch = 16)
lines(10^mod_percentbinding@xCurve, mod_percentbinding@yCurve)
points(absconc_pgmL, pbb0[desesamples], col = "red", pch = 17)

legend(600, 0.75, c("standard curve", "samples"), lty = c(1, NA), col = c("black", "red"), pch = c(16, 16), bty = 'n')
# axis(4, at = pbb0[desesamples], labels = ss$lab[desesamples], las = 1)
text(absconc_pgmL, pbb0[desesamples], ss$lab[desesamples])


orig_pgmL <- absconc_pgmL / ss$dilution[desesamples]
orig_ng_g <- orig_pgmL * (0.250 /1) * (1/1000) * (1/ss$neat_mass_mg[desesamples]) * (1000/1)
# (ng / g) = (pg / mL) * (0.250 mL / 1) * (1 ng / 1000 pg) * (1 / sample mass mg) * (1000 mg / 1 g)

data.frame(ss$lab[desesamples], orig_ng_g, absconc_pgmL)
# write.table(
data.frame(absconc_pgmL, pbb0[desesamples], ss$lab[desesamples])
# ,
	# file = "bemncortplate.csv",
	# row.names = FALSE,
	# sep = ',')