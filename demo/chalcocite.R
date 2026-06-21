# CHNOSZ/demo/chalcocite.R
# Solubility of chalcocite, after Trofimov et al. (2023)
# doi:10.1016/j.chemgeo.2023.121413
# 20260618 first version

library(CHNOSZ)

# Define conditions
T <- 400
P <- 1000
log_mCl <- c(-6, 1.5)

# Assign values for the transect
res <- 40
Cl_vals <- seq(log_mCl[1], log_mCl[2], length.out = res)

# Create the pH buffer
mod.buffer("KMQ", c("quartz", "muscovite", "K-feldspar"), "cr", 0)
# Create the aH2S buffer
mod.buffer("PHM", c("pyrite", "hematite", "magnetite"), "cr", 0)

# Calculate ionic strength from total Cl molality
salt <- lapply(10^Cl_vals, NaCl, T = T, P = P)
# Get it into more convenient format
salt <- as.data.frame(do.call(rbind, lapply(salt, unlist)))

# Calculate pH from KMQ buffer
# Initialize pH
pH <- rep(5, res)
# Iterate until convergence
for(i in 1:5) {

  # Calculate K+ from K-feldspar solubility
  # NOTE 1: Any K-bearing species must be first
  # NOTE 2: muscovite and quartz represent the KMQ buffer
  basis(c("K+", "muscovite", "quartz", "Cl-", "H2O", "oxygen", "H+"))
  # Get aqueous K species
  iaq <- retrieve("K", c("H", "O", "Cl"), state = "aq")
  species("K-feldspar")
  s <- solubility(iaq, "Cl-" = log10(salt$m_Clminus), IS = salt$IS, pH = pH, T = T, P = P)
  log_mK <- s$loga.balance

  # Now calculate pH with the KMQ buffer using the K+ activity we just calculated
  basis(c("K+", "Al2O3", "SiO2", "H2O", "oxygen", "H+"))
  basis("H+", "KMQ")
  # FIXME: using IS here decreases pH and pushes solubilty too high
  #a <- affinity(IS = 10^log_mK, "K+" = log_mK, IS = salt$IS, T = T, P = P, return.buffer = TRUE)
  a <- affinity(IS = 10^log_mK, "K+" = log_mK, T = T, P = P, return.buffer = TRUE)
  # This is our new pH
  pH <- -a$`H+`
  # Print the ranges of K+ activity and pH
  range_log_mK <- paste(round(range(log_mK), 2), collapse = " to ")
  print(paste("log mK+ range:", range_log_mK))
  range_pH <- paste(round(range(pH), 2), collapse = " to ")
  print(paste("pH range:", range_pH))

}

for(substrate in c("chalcopyrite", "gold")) {

  # Get the metal
  metal <- switch(substrate, chalcopyrite = "Cu", gold = "Au")

  # Set up system
  # NOTE 1: Any Cu- or Au- bearing species must be first otherwise we get this cryptic error message from solubility:
  # Error in put.basis(ispecies, logact) :
  #   the number of basis species is greater than the number of elements and charge
  # NOTE 2: Use hematite or magnetite for the Fe component
  #   (to buffer Fe activity, which affects chalcopyrite stability)
  # NOTE 3: Because the activity of H2S is buffered, no mosaic calculation is possible
  basis(c(metal, "magnetite", "Cl-", "H2S", "H2O", "oxygen", "H+"))
  # Set the buffers
  basis("O2", "HM")
  basis("H2S", "PHM")
  ## Don't set KMQ here becuase we'll insert the pH values calculated above
  #basis("H+", "KMQ")

  # Get aqueous metal species
  iaq <- retrieve(metal, c("H", "O", "Cl", "S"), state = "aq")

  # Calculate solubility of substrate
  species(substrate)
  s <- solubility(iaq, "Cl-" = log10(salt$m_Clminus), IS = salt$IS, pH = pH, T = T, P = P)

  # Insert original Cl molality for plotting
  s$vals[[1]] <- Cl_vals

  # Make diagram
  if(substrate == "chalcopyrite") {
    diagram(s, type = "loga.balance", col = 4, lwd = 2, ylim = c(-8, 0))
    diagram(s, type = "loga.equil", add = TRUE)
    ilab <- 25
    text(Cl_vals[ilab], s$loga.balance[ilab] + 0.8, "Cu_tot", font = 2, col = 4)
  } else {
    diagram(s, type = "loga.balance", col = 2, lwd = 2, add = TRUE)
    ilab <- 30
    text(Cl_vals[ilab], s$loga.balance[ilab] + 0.3, "Au_tot", font = 2, col = 2)
    ## Uncomment to plot Au species
    #diagram(s, type = "loga.equil", col = 2, add = TRUE)
  }

}

# Add legend
l_T <- lT(T)
l_P <- lP(P)
l_O2 <- describe.basis()[6]
l_H2S <- describe.basis()[4]
l_pH <- "pH = KMQ"
l_text <- c(l_T, l_P, l_O2, l_H2S, l_pH)
legend("topleft", legend = l_text, bty = "n")
# Add title
title("Chalcocite and gold solubility in rock-buffered system\nAfter Trofimov et al. (2023)", font.main = 1)
