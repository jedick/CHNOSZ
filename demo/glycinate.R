# CHNOSZ/demo/glycinate.R
# Plot logK of metal-glycinate complexes
# 20190207

library(CHNOSZ)

# Divalent metals
di <- c("Cu+2", "Ni+2", "Co+2", "Mn+2", "Zn+2", "Cd+2")
# Divalent metals with one glycinate
di1 <- c("Cu(Gly)+", "Ni(Gly)+", "Co(Gly)+", "Mn(Gly)+", "Zn(Gly)+", "Cd(Gly)+")
# Divalent metals with two glycinates
di2 <- c("Cu(Gly)2", "Ni(Gly)2", "Co(Gly)2", "Mn(Gly)2", "Zn(Gly)2", "Cd(Gly)2")
# Monovalent metals
mo <- c("Au+", "Ag+", "Na+", "Tl+", "Cu+")
# Monovalent metals with one glycinate
mo1 <- c("Au(Gly)", "Ag(Gly)", "Na(Gly)", "Tl(Gly)", "Cu(Gly)")
# Monovalent metals with two glycinates
mo2 <- c("Au(Gly)2-", "Ag(Gly)2-", "Na(Gly)2-", "Tl(Gly)2-", "Cu(Gly)2-")

# Set the temperature values
T <- seq(0, 150, 10)
# Calculate the logKs using data from Azadi et al., 2019
# doi:10.1016/j.fluid.2018.10.002
logK_di1 <- logK_di2 <- logK_mo1 <- logK_mo2 <- list()
for(i in 1:length(di1)) logK_di1[[i]] <- subcrt(c(di[i], "glycinate", di1[i]), c(-1, -1, 1), T = T)$out$logK
for(i in 1:length(di2)) logK_di2[[i]] <- subcrt(c(di[i], "glycinate", di2[i]), c(-1, -2, 1), T = T)$out$logK
for(i in 1:length(mo1)) logK_mo1[[i]] <- subcrt(c(mo[i], "glycinate", mo1[i]), c(-1, -1, 1), T = T)$out$logK
for(i in 1:length(mo2)) logK_mo2[[i]] <- subcrt(c(mo[i], "glycinate", mo2[i]), c(-1, -2, 1), T = T)$out$logK

# Calculate the logKs for divalent metals using data from Shock and Koretsky, 1995
# doi:10.1016/0016-7037(95)00058-8
add.OBIGT(system.file("extdata/misc/SK95.csv", package = "CHNOSZ"))
logK_di1_SK95 <- logK_di2_SK95 <- list()
for(i in 1:length(di1)) logK_di1_SK95[[i]] <- subcrt(c(di[i], "glycinate", di1[i]), c(-1, -1, 1), T = T)$out$logK
for(i in 1:length(di2)) logK_di2_SK95[[i]] <- subcrt(c(di[i], "glycinate", di2[i]), c(-1, -2, 1), T = T)$out$logK
reset()

# Set up the plots
opar <- par(no.readonly = TRUE)
layout(matrix(1:6, byrow = TRUE, nrow = 2), widths = c(2, 2, 1))
par(mar = c(4, 3.2, 2.5, 0.5), mgp = c(2.1, 1, 0), las = 1, cex = 0.8)
xlab <- axis.label("T")
ylab <- axis.label("logK")

# First row: divalent metals
matplot(T, sapply(logK_di1, c), type = "l", lwd = 2, lty = 1, xlab = xlab, ylab = ylab)
matplot(T, sapply(logK_di1_SK95, c), type = "l", lwd = 2, lty = 2, add = TRUE)
legend(-9, 7.7, c("Azadi et al., 2019", "Shock and Koretsky, 1995"), lty = c(1, 2), bty = "n", cex = 1)
mtext(expression(M^"+2" + Gly^"-" == M*(Gly)^"+"), line = 0.5)
matplot(T, sapply(logK_di2, c), type = "l", lwd = 2, lty = 1, xlab = xlab, ylab = ylab)
matplot(T, sapply(logK_di2_SK95, c), type = "l", lwd = 2, lty = 2, add = TRUE)
legend(-9, 14, c("Azadi et al., 2019", "Shock and Koretsky, 1995"), lty = c(1, 2), bty = "n", cex = 1)
mtext(expression(M^"+2" + 2*Gly^"-" == M*(Gly)[2]), line = 0.5)
plot.new()
par(xpd = NA)
legend("right", as.expression(lapply(di, expr.species)), lty = 1, col = 1:6, bty = "n", cex = 1.2, lwd = 2)

# Add overall title
text(0, 1, "metal-\nglycinate\ncomplexes", cex = 1.3, font = 2)
par(xpd = FALSE)

# Second row: monovalent metals
matplot(T, sapply(logK_mo1, c), type = "l", lwd = 2, lty = 1, xlab = xlab, ylab = ylab)
mtext(expression(M^"+" + Gly^"-" == M*(Gly)), line = 0.5)
matplot(T, sapply(logK_mo2, c), type = "l", lwd = 2, lty = 1, xlab = xlab, ylab = ylab)
mtext(expression(M^"+" + 2*Gly^"-" == M*(Gly)[2]^"-"), line = 0.5)
plot.new()
par(xpd = NA)
legend("right", as.expression(lapply(mo, expr.species)), lty = 1, col = 1:5, bty = "n", cex = 1.2, lwd = 2)
par(xpd = FALSE)

layout(matrix(1))
par(opar)
