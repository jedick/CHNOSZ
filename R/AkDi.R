# CHNOSZ/AkDi.R
# Akinfiev-Diamond model for aqueous species
# 20190219 first version

AkDi <- function(property = NULL, parameters = NULL, T = 298.15, P = 1, isPsat = TRUE) {

  # some constants (from Akinfiev and Diamond, 2004 doi:10.1016/j.fluid.2004.06.010)
  MW <- 18.0153 # g mol-1
  NW <- 1000/MW # mol kg-1
  R <- 8.31441 # J K-1 mol-1

  # a list for the output
  out <- list()
  # loop over species
  nspecies <- nrow(parameters)
  for(i in seq_len(nspecies)) {
    PAR <- parameters[i, ]
    # start with an NA-filled data frame
    myprops <- as.data.frame(matrix(NA, ncol=length(property), nrow=length(T)))
    colnames(myprops) <- property
    # just calculate G for now
    for(j in seq_along(property)) {
      if(property[[j]]=="G") {
        # send a message
        message("AkDi(): Akinfiev-Diamond model for ", PAR$name, " gas to aq")
        # get gas properties (J mol-1)
        G_gas <- subcrt(PAR$name, "gas", T=T, P=P, convert=FALSE)$out[[1]]$G
        # TODO: does this work if E.units is cal or J?
        G_gas <- convert(G_gas, "J", T=T)
        # get H2O fugacity (bar)
        GH2O_P <- water("G", T=T, P=P)$G
        GH2O_1 <- water("G", T=T, P=1)$G
        f1 <- exp ( (GH2O_P - GH2O_1) / (1.9872 * T) )
        # for Psat, calculate the real liquid-vapor curve (not 1 bar below 100 degC)
        if(isPsat) {
          P <- water("Psat", T = T, P = "Psat", P1 = FALSE)$Psat
          f1[P < 1] <- P[P < 1]
        }
        # density (g cm-3)
        rho1 <- water("rho", T=T, P=P)$rho / 1000
        # calculate G_hyd (J mol-1)
        G_hyd <- R*T * ( -log(NW) + (1 - PAR$xi) * log(f1) + PAR$xi * log(10 * R * T * rho1 / MW) + rho1 * (PAR$a + PAR$b * (1000/T)^0.5) )
        # calculate the chemical potential (J mol-1)
        G <- G_gas + G_hyd
        # convert J to cal
        G <- convert(G, "cal", T=T)
        # insert into data frame of properties
        myprops$G <- G
      }
    }
    out[[i]] <- myprops
  }
  out
}
