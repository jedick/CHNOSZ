# solubility.R: vectorized solubility calculations without uniroot
# 20181031 jmd

solubility <- function(eout, exp = 1) {
  # exp = 1: e.g. dissolution of CO2
  # exp = 2: e.g. dissolution (dissociation) of CaCO3

  # bookkeeping: track any single species
  itrack <- 1
  # the log activity used to calculate the affinity
  loga.species.track <- eout$species$logact[itrack]
  # the affinities at the starting loga.balance
  A.track <- eout$values[[itrack]]
  # the loga.equil at the starting loga.balance
  loga.equil.track <- eout$loga.equil[[itrack]]

  # subjunctive: what would the affinities be if the
  # activity of the tracked species was set to loga.equil?
  A.whatif <- loga.species.track + A.track - loga.equil.track

  # predictive: assuming the species distribution doesn't change,
  # what is the total loga that would give zero affinity?
  # TODO: modify this according to stoichiometry (species with > 1 of the balanced basis species)
  loga.total <- (eout$loga.balance + A.whatif) / exp

  message("solubility: calculated logarithm of total activity of ", eout$balance)

  # use the predicted loga.total to re-calculate activities of species
  aout <- eout[1:which(names(eout)=="values")]
  equilibrate(aout, loga.balance = loga.total)
}

