# CHNOSZ/demo/solubility.R: vectorized solubility without uniroot
# adapted from CHNOSZ/demo/oldsolub.R
# 20181030 jmd

# for comparison with published calcite solubility plot, see Fig. 4A in
# Manning et al., 2013, Reviews in Mineralogy & Geochemistry, v. 75, pp. 109-148
# (doi: 10.2138/rmg.2013.75.5)

# for comparison with published CO2 solubility plot, see Fig. 4.5 in
# Stumm and Morgan, 1996, Aquatic Chemistry: Chemical Equilibria and Rates in Natural Waters
# (New York: John Wiley & Sons), 3rd edition

par(mfrow=c(1, 2))

for(what in c("CO2", "calcite")) {

  # set up system
  if(what=="CO2") {
    basis("CHNOS+")
    basis("CO2", "gas")
    # ca. atmospheric PCO2
    basis("CO2", -3.5)
  } else if(what=="calcite") {
    basis(c("calcite", "Ca+2", "H2O", "O2", "H+"))
  }
  species(c("CO2", "HCO3-", "CO3-2"))

  # set pH range and resolution, constant temperature and ionic strength
  pH <- c(0, 14)
  res <- 100
  T <- 25
  IS <- 0

  # start with loga.balance = 0
  loga.balance <- 0
  a0 <- affinity(pH = c(pH, res), T = T, IS = IS)
  e0 <- equilibrate(a0, loga.balance = loga.balance)

  # bookkeeping: track any single species
  itrack <- 1
  # the log activity used to calculate the affinity
  loga.species.track <- species(itrack)$logact
  # the affinities at the starting loga.balance
  A.track <- a0$values[[itrack]]
  # the loga.equil at the starting loga.balance
  loga.equil.track <- e0$loga.equil[[itrack]]

  # subjunctive: what would the affinities be if the
  # activity of the tracked species was set to loga.equil?
  A.whatif <- loga.species.track + A.track - loga.equil.track

  # predictive: assuming the species distribution doesn't change,
  # what is the log(total activity) that gives zero affinity?
  # TODO: modify this according to stoichiometry (species with > 1 of the balanced basis species)
  loga.total <- loga.balance + A.whatif

  # a dissociation reaction makes two things, so the exponent is not 1
  if(what=="calcite") loga.total <- loga.total / 2

  # use the predicted loga.total to re-calculate activities of species
  e1 <- equilibrate(a0, loga.balance = loga.total)

#  # check that we got stable conditions
#  # (not easily vectorized because we change the activities of species)
#  print("checking for stable conditions (affinity = 0)")
#  for(i in 1:length(a0$vals[[1]])) {
#    basis(a0$vars[1], a0$vals[[1]][i])
#    if(what=="calcite") basis("Ca+2", loga.total[i])
#    species(1:3, sapply(e1$loga.equil, "[", i))
#    atest <- suppressMessages(affinity(T = T, IS = IS))
#    stopifnot(all(sapply(as.numeric(unlist(atest$values)), all.equal, 0)))
#  }

  # make plot
  ylim <- c(-10, 4)
  thermo.plot.new(xlim = range(pH), ylim = ylim, xlab = "pH", ylab = "log a")
  lines(a0$vals[[1]], loga.total, lwd = 4, col = "green2")
  lines(a0$vals[[1]], e1$loga.equil[[1]], lwd = 2)
  lines(a0$vals[[1]], e1$loga.equil[[2]], lty = 2, lwd = 2)
  lines(a0$vals[[1]], e1$loga.equil[[3]], lty = 3, lwd = 2)

  legend(ifelse(what=="calcite", "topright", "topleft"), lty = c(1, 1:3), lwd = c(4, 2, 2, 2), col = c("green2", rep("black", 3)),
         legend = as.expression(c("total", expr.species("CO2", state = "aq"), expr.species("HCO3-"), expr.species("CO3-2"))))
  title(main = substitute("Solubility of"~what~"at"~T~degree*"C", list(what = what, T = T)))

}
