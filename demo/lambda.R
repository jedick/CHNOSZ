# plot effects of lambda transition in quartz
# after Berman 1988 Figs. 1 and 2
library(CHNOSZ)

opar <- par(no.readonly = TRUE)
layout(matrix(c(1, 4:2, 1, 7:5), nrow=4), heights=c(0.7, 3, 3, 3))
# plot title first
par(mar=c(0, 0, 0, 0))
plot.new()
text(0.5, 0.5, "Effects of lambda transition in quartz, after Berman (1988) Figs. 1 and 2", cex=1.8)
par(mar=c(4, 4.5, 1, 0.5), cex=0.8)

TC <- 0:1200
T <- convert(TC, "K")

labplot <- function(x) label.plot(x, xfrac=0.9, yfrac=0.1, paren=TRUE)
Cplab <- axis.label("Cp")
Vlab <- axis.label("V")
Tlab <- axis.label("T")

# calculate properties at 1 kbar with and without transition
Qz_1bar <- Berman("quartz", T=T, units="J")
Qz_1bar_notrans <- Berman("quartz", T=T, calc.transition=FALSE, units="J")

## Fig. 1a: volume
plot(TC, Qz_1bar$V, type="l", xlab=Tlab, ylab=Vlab, ylim=c(22.5, 24))
# FIXME: why don't we get the curvature his plot for V shows?
# Should it be in the v4 parameter (but it's zero)??
## add data points digitized from Fig. 3B of Helgeson et al., 1978
## and Fig. 1a of Berman, 1988
skinner <- list(T=c(550, 560, 570, 580, 590), V=c(23.44, 23.48, 23.54, 23.72, 23.72))
robie <- list(T=575, V=23.72)
ghioroso <- list(T=c(0, 100, 200, 300, 400, 500, 575, 575, 675, 775, 875, 975),
                 V=c(22.7, 22.77, 22.85, 22.94, 23.1, 23.3, 23.5, 23.71, 23.69, 23.69, 23.71, 23.73))
points(skinner$T, skinner$V, pch=19)
points(robie$T, robie$V)
points(ghioroso$T, ghioroso$V, pch=0)

## calculate finite difference derivative of dGlambda/dP = Vlambda between 1 and 1.001 bar
Glambda_1bar <- Qz_1bar_notrans$G - Qz_1bar$G
Qz_1.001bar <- Berman("quartz", T=T, P=1.001, units="J")
Qz_1.001bar_notrans <- Berman("quartz", T=T, P=1.001, calc.transition=FALSE, units="J")
Glambda_1.001bar <- Qz_1.001bar_notrans$G - Qz_1.001bar$G
dGlambdadP <- (Glambda_1.001bar - Glambda_1bar) / 0.001
# we're using Joules, so multiply by ten to get cm^3
Vlambda <- -10 * dGlambdadP 
VQz <- Qz_1bar$V + Vlambda
# above 848K, we have beta quartz
Qz_beta <- Berman("quartz,beta", T=T, P=1, units="J")
VQz[T >= 848.15] <- Qz_beta$V[T >= 848.14]
lines(TC, VQz, lty=2)
legend("topleft", legend="1 bar", bty="n")
labplot("a")

## Fig. 1b: heat capacity
plot(TC, Qz_1bar$Cp, type="l", xlab=Tlab, ylab=Cplab)
lines(TC, Qz_1bar_notrans$Cp, lty=3)
legend("topleft", legend="1 bar", bty="n")
labplot("b")

## calculate properties at 10 kbar with and without transition
Qz_10bar <- Berman("quartz", T=T, P=10000, units="J")
Qz_10bar_notrans <- Berman("quartz", T=T, P=10000, calc.transition=FALSE, units="J")
## Fig. 1c: heat capacity
plot(TC, Qz_10bar$Cp, type="l", xlab=Tlab, ylab=Cplab)
lines(TC, Qz_10bar_notrans$Cp, lty=3)
legend("topleft", legend="10 kb", bty="n")
labplot("c")

# like Ber88 Fig. 2
Tlambda <- 848 # Kelvin
dTdP <- 0.0237
Pkb <- seq(1, 50, 1)
P <- 1000 * Pkb
T <- Tlambda + dTdP * (P - 1)
Qz_withtrans <- Berman("quartz", T=T, P=P, units="J")
Qz_notrans <- Berman("quartz", T=T, P=P, calc.transition=FALSE, units="J")
Qz_lambda <- Qz_withtrans - Qz_notrans
Plab <- expression(list(italic(P), "kb"))
plot(Pkb, Qz_lambda$G, type="l", ylim=c(-300, -50), ylab=axis.label("DlG"), xlab=Plab)
labplot("d")
plot(Pkb, Qz_lambda$H, type="l", ylim=c(1200, 1800), ylab=axis.label("DlH"), xlab=Plab)
labplot("e")
plot(Pkb, Qz_lambda$S, type="l", ylim=c(0, 3), ylab=axis.label("DlS"), xlab=Plab)
labplot("f")

reset()

layout(matrix(1))
par(opar)
