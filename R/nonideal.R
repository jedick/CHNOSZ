# CHNOSZ/nonideal.R
# first version of function: 20080308 jmd
# moved to nonideal.R from util.misc.R 20151107
# added Helgeson method 20171012

nonideal <- function(species, speciesprops, IS, T, P, A_DH, B_DH, m_star=NULL, method=thermo()$opt$nonideal) {
  # Generate nonideal contributions to thermodynamic properties
  # number of species, same length as speciesprops list
  # T in Kelvin, same length as nrows of speciespropss
  # arguments A_DH and B_DH are needed for all methods other than "Alberty", and P is needed for "bgamma"
  # m_start is the total molality of all dissolved species; if not given, it is taken to be equal to ionic strength

  mettext <- function(method) {
    mettext <- paste(method, "equation")
    if(method=="Bdot0") mettext <- "B-dot equation (B-dot = 0)"
    mettext
  }

  # We can use this function to change the nonideal method option
  if(missing(speciesprops)) {
    if(species[1] %in% c("Bdot", "Bdot0", "bgamma", "bgamma0", "Alberty")) {
      thermo <- get("thermo", CHNOSZ)
      oldnon <- thermo$opt$nonideal
      thermo$opt$nonideal <- species[1]
      assign("thermo", thermo, CHNOSZ)
      message("nonideal: setting nonideal option to use ", mettext(species))
      return(invisible(oldnon))
    } else stop(species[1], " is not a valid nonideality setting (Bdot, Bdot0, bgamma, bgamma0, or Alberty)")
  }

  # Check if we have a valid method setting
  if(!method %in% c("Alberty", "Bdot", "Bdot0", "bgamma", "bgamma0")) {
    if(missing(method)) stop("invalid setting (", thermo$opt$nonideal, ") in thermo()$opt$nonideal")
    else stop("invalid method (", thermo$opt$nonideal, ")")
  }

  #R <- 1.9872  # gas constant, cal K^-1 mol^-1
  R <- 8.314445  # gas constant, J K^-1 mol^-1  20220325

  # Function to calculate extended Debye-Huckel equation and derivatives using Alberty's parameters
  Alberty <- function(prop = "loggamma", Z, I, T) {
    # Extended Debye-Huckel equation ("log")
    # and its partial derivatives ("G","H","S","Cp")
    # T in Kelvin
    B <- 1.6 # L^0.5 mol^-0.5 (Alberty, 2003 p. 47)
    # Equation for A from Clarke and Glew, 1980
    #A <- expression(-16.39023 + 261.3371/T + 3.3689633*log(T)- 1.437167*(T/100) + 0.111995*(T/100)^2)
    # A = alpha / 3 (Alberty, 2001)
    alpha <- expression(3 * (-16.39023 + 261.3371/T + 3.3689633*log(T)- 1.437167*(T/100) + 0.111995*(T/100)^2))
    ## Equation for alpha from Alberty, 2003 p. 48
    #alpha <- expression(1.10708 - 1.54508E-3 * T + 5.95584E-6 * T^2)
    # from examples for deriv() to take first and higher-order derivatives
    DD <- function(expr, name, order = 1) {
      if(order < 1) stop("'order' must be >= 1")
      if(order == 1) D(expr, name)
      else DD(D(expr, name), name, order - 1)
    }
    # Alberty, 2003 Eq. 3.6-1
    lngamma <- function(alpha, Z, I, B) - alpha * Z^2 * I^(1/2) / (1 + B * I^(1/2))
    # 20171013 convert lngamma to common logarithm
    # 20190603 use equations for H, S, and Cp from Alberty, 2001 (doi:10.1021/jp011308v)
    if(prop=="loggamma") return(lngamma(eval(alpha), Z, I, B) / log(10))
    else if(prop=="G") return(R * T * lngamma(eval(alpha), Z, I, B))
    else if(prop=="H") return(- R * T^2 * lngamma(eval(DD(alpha, "T", 1)), Z, I, B))
    else if(prop=="S") return( ( - R * T^2 * lngamma(eval(DD(alpha, "T", 1)), Z, I, B) - R * T * lngamma(eval(alpha), Z, I, B) ) / T)
    else if(prop=="Cp") return(- 2 * R * T * lngamma(eval(DD(alpha, "T", 1)), Z, I, B) - R * T^2 * lngamma(eval(DD(alpha, "T", 2)), Z, I, B))
  }
  
  # Function for Debye-Huckel equation with b_gamma or B-dot extended term parameter (Helgeson, 1969)
  Helgeson <- function(prop = "loggamma", Z, I, T, A_DH, B_DH, acirc, m_star, bgamma) {
    loggamma <- - A_DH * Z^2 * I^0.5 / (1 + acirc * B_DH * I^0.5) - log10(1 + 0.0180153 * m_star) + bgamma * I
    if(prop=="loggamma") return(loggamma)
    else if(prop=="G") return(R * T * log(10) * loggamma)
    # note the log(10) (=2.303) ... use natural logarithm to calculate G
  }

  # Function for Setchenow equation with b_gamma or B-dot extended term parameter (Shvarov and Bastrakov, 1999)  20181106
  Setchenow <- function(prop = "loggamma", I, T, m_star, bgamma) {
    loggamma <- - log10(1 + 0.0180153 * m_star) + bgamma * I
    if(prop=="loggamma") return(loggamma)
    else if(prop=="G") return(R * T * log(10) * loggamma)
  }

  # Get species indices
  if(!is.numeric(species[[1]])) species <- info(species, "aq")
  # loop over species #1: get the charge
  Z <- numeric(length(species))
  for(i in 1:length(species)) {
    # force a charge count even if it's zero
    mkp <- makeup(c("Z0", species[i]), sum=TRUE)
    thisZ <- mkp[match("Z", names(mkp))]
    # no charge if Z is absent from the formula or equal to zero
    if(is.na(thisZ)) next
    if(thisZ==0) next
    Z[i] <- thisZ
  }
  # get species formulas to assign acirc 20181105
  formula <- get("thermo", CHNOSZ)$OBIGT$formula[species]
  if(grepl("Bdot", method)) {
    # "ion size paramter" taken from UT_SIZES.REF of HCh package (Shvarov and Bastrakov, 1999),
    # based on Table 2.7 of Garrels and Christ, 1965
    acircdat <- c("Rb+"=2.5, "Cs+"=2.5, "NH4+"=2.5, "Tl+"=2.5, "Ag+"=2.5,
      "K+"=3, "Cl-"=3, "Br-"=3, "I-"=3, "NO3-"=3,
      "OH-"=3.5, "F-"=3.5, "HS-"=3.5, "BrO3-"=3.5, "IO3-"=3.5, "MnO4-"=3.5,
      "Na+"=4, "HCO3-"=4, "H2PO4-"=4, "HSO3-"=4, "Hg2+2"=4, "SO4-2"=4, "SeO4-2"=4, "CrO4-2"=4, "HPO4-2"=4, "PO4-3"=4,
      "Pb+2"=4.5, "CO3-2"=4.5, "SO4-2"=4.5, "MoO4-2"=4.5,
      "Sr+2"=5, "Ba+2"=5, "Ra+2"=5, "Cd+2"=5, "Hg+2"=5, "S-2"=5, "WO4-2"=5,
      "Li+"=6, "Ca+2"=6, "Cu+2"=6, "Zn+2"=6, "Sn+2"=6, "Mn+2"=6, "Fe+2"=6, "Ni+2"=6, "Co+2"=6,
      "Mg+2"=8, "Be+2"=8,
      "H+"=9, "Al+3"=9, "Cr+3"=9, "La+3"=9, "Ce+3"=9, "Y+3"=9, "Eu+3"=9,
      "Th+4"=11, "Zr+4"=11, "Ce+4"=11, "Sn+4"=11)
    acirc <- as.numeric(acircdat[formula])
    acirc[is.na(acirc)] <- 4.5
    ## Make a message
    #nZ <- sum(Z!=0)
    #if(nZ > 1) message("nonideal: using ", paste(acirc[Z!=0], collapse=" "), " for ion size parameters of ", paste(formula[Z!=0], collapse=" "))
    #else if(nZ==1) message("nonideal: using ", acirc[Z!=0], " for ion size parameter of ", formula[Z!=0])
    # use correct units (cm) for ion size parameter
    acirc <- acirc * 10^-8
  } else if(grepl("bgamma", method)) {
    # "distance of closest approach" of ions in NaCl solutions (HKF81 Table 2)
    acirc <- rep(3.72e-8, length(species))
  }
  # Get b_gamma or B-dot
  if(method=="bgamma") bgamma <- bgamma(convert(T, "C"), P)
  else if(method=="Bdot") bgamma <- Bdot(convert(T, "C"))
  else if(method %in% c("Bdot0", "bgamma0")) bgamma <- 0
  # Loop over species #2: activity coefficient calculations
  if(is.null(m_star)) m_star <- IS
  iH <- info("H+")
  ie <- info("e-")
  speciesprops <- as.list(speciesprops)
  icharged <- ineutral <- logical(length(species))
  for(i in 1:length(species)) {
    myprops <- speciesprops[[i]]
    # To keep unit activity coefficients of the proton and electron
    if(species[i] == iH & get("thermo", CHNOSZ)$opt$ideal.H) next
    if(species[i] == ie & get("thermo", CHNOSZ)$opt$ideal.e) next
    didcharged <- didneutral <- FALSE
    # Logic for neutral and charged species 20181106
    if(Z[i]==0) {
      for(j in 1:ncol(myprops)) {
        pname <- colnames(myprops)[j]
        if(!pname %in% c("G", "H", "S", "Cp")) next
        if(identical(get("thermo", CHNOSZ)$opt$Setchenow, "bgamma")) {
          myprops[, j] <- myprops[, j] + Setchenow(pname, IS, T, m_star, bgamma)
          didneutral <- TRUE
        } else if(identical(get("thermo", CHNOSZ)$opt$Setchenow, "bgamma0")) {
          myprops[, j] <- myprops[, j] + Setchenow(pname, IS, T, m_star, bgamma = 0)
          didneutral <- TRUE
        }
      }
    } else {
      for(j in 1:ncol(myprops)) {
        pname <- colnames(myprops)[j]
        if(!pname %in% c("G", "H", "S", "Cp")) next
        if(method=="Alberty") {
          myprops[, j] <- myprops[, j] + Alberty(pname, Z[i], IS, T)
          didcharged <- TRUE
        } else {
          myprops[, j] <- myprops[, j] + Helgeson(pname, Z[i], IS, T, A_DH, B_DH, acirc[i], m_star, bgamma)
          didcharged <- TRUE
        }
      }
    }
    # Append a loggam column if we did any calculations of adjusted thermodynamic properties
    if(didcharged) {
      if(method=="Alberty") myprops <- cbind(myprops, loggam = Alberty("loggamma", Z[i], IS, T))
      else myprops <- cbind(myprops, loggam = Helgeson("loggamma", Z[i], IS, T, A_DH, B_DH, acirc[i], m_star, bgamma))
    }
    if(didneutral) {
      if(get("thermo", CHNOSZ)$opt$Setchenow == "bgamma") myprops <- cbind(myprops, loggam = Setchenow("loggamma", IS, T, m_star, bgamma))
      else if(get("thermo", CHNOSZ)$opt$Setchenow == "bgamma0") myprops <- cbind(myprops, loggam = Setchenow("loggamma", IS, T, m_star, bgamma = 0))
    }
    # Save the calculated properties and increment progress counters
    speciesprops[[i]] <- myprops
    if(didcharged) icharged[i] <- TRUE
    if(didneutral) ineutral[i] <- TRUE
  }
  if(sum(icharged) > 0) message("nonideal: calculations for ", paste(formula[icharged], collapse=", "), " (", mettext(method), ")")
  if(sum(ineutral) > 0) message("nonideal: calculations for ", paste(formula[ineutral], collapse=", "), " (Setchenow equation)")
  return(speciesprops)
}

bgamma <- function(TC = 25, P = 1, showsplines = "") {
  # 20171012 calculate b_gamma using P, T, points from:
  # Helgeson, 1969 (doi:10.2475/ajs.267.7.729)
  # Helgeson et al., 1981 (doi:10.2475/ajs.281.10.1249)
  # Manning et al., 2013 (doi:10.2138/rmg.2013.75.5)
  # T in degrees C
  T <- TC
  # are we at a pre-fitted constant pressure?
  uP <- unique(P)
  is1 <- identical(uP, 1) & all(T==25)
  is500 <- identical(uP, 500)
  is1000 <- identical(uP, 1000)
  is2000 <- identical(uP, 2000)
  is3000 <- identical(uP, 3000)
  is4000 <- identical(uP, 4000)
  is5000 <- identical(uP, 5000)
  is10000 <- identical(uP, 10000)
  is20000 <- identical(uP, 20000)
  is30000 <- identical(uP, 30000)
  is40000 <- identical(uP, 40000)
  is50000 <- identical(uP, 50000)
  is60000 <- identical(uP, 60000)
  isoP <- is1 | is500 | is1000 | is2000 | is3000 | is4000 | is5000 | is10000 | is20000 | is30000 | is40000 | is50000 | is60000
  # values for Bdot x 100 from Helgeson (1969), Figure (P = Psat)
  if(!isoP | showsplines != "") {
    T0 <- c(23.8, 49.4, 98.9, 147.6, 172.6, 197.1, 222.7, 248.1, 268.7)
    B0 <- c(4.07, 4.27, 4.30, 4.62, 4.86, 4.73, 4.09, 3.61, 1.56) / 100
    # we could use the values from Hel69 Table 2 but extrapolation of the
    # their fitted spline function turns sharply upward above 300 degC
    #T0a <- c(25, 50, 100, 150, 200, 250, 270, 300)
    #B0a <- c(4.1, 4.35, 4.6, 4.75, 4.7, 3.4, 1.5, 0)
    S0 <- splinefun(T0, B0)
  }
  # values for bgamma x 100 from Helgeson et al., 1981 Table 27 
  if(is500 | !isoP | showsplines != "") {
    T0.5 <- seq(0, 400, 25)
    B0.5 <- c(5.6, 7.1, 7.8, 8.0, 7.8, 7.5, 7.0, 6.4, 5.7, 4.8, 3.8, 2.6, 1.0, -1.2, -4.1, -8.4, -15.2) / 100
    S0.5 <- splinefun(T0.5, B0.5)
    if(is500) return(S0.5(T))
  }
  if(is1000 | !isoP | showsplines != "") {
    T1 <- seq(0, 500, 25)
    B1 <- c(6.6, 7.7, 8.7, 8.3, 8.2, 7.9, 7.5, 7.0, 6.5, 5.9, 5.2, 4.4, 3.5, 2.5, 1.1, -0.6, -2.8, -5.7, -9.3, -13.7, -19.2) / 100
    S1 <- splinefun(T1, B1)
    if(is1000) return(S1(T))
  }
  if(is2000 | !isoP | showsplines != "") {
    # 550 and 600 degC points from Manning et al., 2013 Fig. 11
    T2 <- c(seq(0, 500, 25), 550, 600)
    B2 <- c(7.4, 8.3, 8.8, 8.9, 8.9, 8.7, 8.5, 8.1, 7.8, 7.4, 7.0, 6.6, 6.2, 5.8, 5.2, 4.6, 3.8, 2.9, 1.8, 0.5, -1.0, -3.93, -4.87) / 100
    S2 <- splinefun(T2, B2)
    if(is2000) return(S2(T))
  }
  if(is3000 | !isoP | showsplines != "") {
    T3 <- seq(0, 500, 25)
    B3 <- c(6.5, 8.3, 9.2, 9.6, 9.7, 9.6, 9.4, 9.3, 9.2, 9.0, 8.8, 8.6, 8.3, 8.1, 7.8, 7.5, 7.1, 6.6, 6.0, 5.4, 4.8) / 100
    S3 <- splinefun(T3, B3)
    if(is3000) return(S3(T))
  }
  if(is4000 | !isoP | showsplines != "") {
    T4 <- seq(0, 500, 25)
    B4 <- c(4.0, 7.7, 9.5, 10.3, 10.7, 10.8, 10.8, 10.8, 10.7, 10.6, 10.5, 10.4, 10.3, 10.2, 10.0, 9.8, 9.6, 9.3, 8.9, 8.5, 8.2) / 100
    S4 <- splinefun(T4, B4)
    if(is4000) return(S4(T))
  }
  if(is5000 | !isoP | showsplines != "") {
    # 550 and 600 degC points from Manning et al., 2013 Fig. 11
    T5 <- c(seq(0, 500, 25), 550, 600)
    B5 <- c(0.1, 6.7, 9.6, 11.1, 11.8, 12.2, 12.4, 12.4, 12.4, 12.4, 12.4, 12.3, 12.3, 12.2, 12.1, 11.9, 11.8, 11.5, 11.3, 11.0, 10.8, 11.2, 12.52) / 100
    S5 <- splinefun(T5, B5)
    if(is5000) return(S5(T))
  }
  # 10, 20, and 30 kb points from Manning et al., 2013 Fig. 11
  # here, one control point at 10 degC is added to make the splines curve down at low T
  if(is10000 | !isoP | showsplines != "") {
    T10 <- c(25, seq(300, 1000, 50))
    B10 <- c(12, 17.6, 17.8, 18, 18.2, 18.9, 21, 23.3, 26.5, 28.8, 31.4, 34.1, 36.5, 39.2, 41.6, 44.1) / 100
    S10 <- splinefun(T10, B10)
    if(is10000) return(S10(T))
  }
  if(is20000 | !isoP | showsplines != "") {
    T20 <- c(25, seq(300, 1000, 50))
    B20 <- c(16, 21.2, 21.4, 22, 22.4, 23.5, 26.5, 29.2, 32.6, 35.2, 38.2, 41.4, 44.7, 47.7, 50.5, 53.7) / 100
    S20 <- splinefun(T20, B20)
    if(is20000) return(S20(T))
  }
  if(is30000 | !isoP | showsplines != "") {
    T30 <- c(25, seq(300, 1000, 50))
    B30 <- c(19, 23.9, 24.1, 24.6, 25.2, 26.7, 30.3, 32.9, 36.5, 39.9, 43, 46.4, 49.8, 53.2, 56.8, 60) / 100
    S30 <- splinefun(T30, B30)
    if(is30000) return(S30(T))
  }
  # 40-60 kb points extrapolated from 10-30 kb points of Manning et al., 2013
  if(is40000 | !isoP | showsplines != "") {
    T40 <- c(seq(300, 1000, 50))
    B40 <- c(25.8, 26, 26.4, 27.2, 28.9, 33, 35.5, 39.2, 43.2, 46.4, 49.9, 53.4, 57.1, 61.2, 64.4) / 100
    S40 <- splinefun(T40, B40)
    if(is40000) return(S40(T))
  }
  if(is50000 | !isoP | showsplines != "") {
    T50 <- c(seq(300, 1000, 50))
    B50 <- c(27.1, 27.3, 27.7, 28.5, 30.5, 34.8, 37.3, 41.1, 45.5, 48.7, 52.4, 55.9, 59.8, 64.3, 67.5) / 100
    S50 <- splinefun(T50, B50)
    if(is50000) return(S50(T))
  }
  if(is60000 | !isoP | showsplines != "") {
    T60 <- c(seq(300, 1000, 50))
    B60 <- c(28, 28.2, 28.6, 29.5, 31.6, 36.1, 38.6, 42.5, 47.1, 50.4, 54.1, 57.6, 61.6, 66.5, 69.7) / 100
    S60 <- splinefun(T60, B60)
    if(is60000) return(S60(T))
  }
  # show points and spline(T) curves
  if(showsplines == "T") {
    thermo.plot.new(c(0, 1000), c(-.2, .7), xlab=axis.label("T"), ylab=expression(italic(b)[gamma]))
    points(T0, B0, pch=0)
    points(T0.5, B0.5, pch=1)
    points(T1, B1, pch=1)
    points(T2[-c(22:23)], B2[-c(22:23)], pch=1)
    points(T2[c(22:23)], B2[c(22:23)], pch=2)
    points(T3, B3, pch=1)
    points(T4, B4, pch=1)
    points(T5[-c(22:23)], B5[-c(22:23)], pch=1)
    points(T5[c(22:23)], B5[c(22:23)], pch=2)
    points(T10[-1], B10[-1], pch=2)
    points(T20[-1], B20[-1], pch=2)
    points(T30[-1], B30[-1], pch=2)
    points(T10[1], B10[1], pch=5)
    points(T20[1], B20[1], pch=5)
    points(T30[1], B30[1], pch=5)
    points(T40, B40, pch=6)
    points(T50, B50, pch=6)
    points(T60, B60, pch=6)
    col <- rev(topo.colors(13))
    T0 <- seq(0, 350, 5); lines(T0, S0(T0), col=col[1])
    T0.5 <- seq(0, 500, 5); lines(T0.5, S0.5(T0.5), col=col[2])
    T1 <- seq(0, 500, 5); lines(T1, S1(T1), col=col[3])
    T2 <- seq(0, 600, 5); lines(T2, S2(T2), col=col[4])
    T3 <- seq(0, 600, 5); lines(T3, S3(T3), col=col[5])
    T4 <- seq(0, 600, 5); lines(T4, S4(T4), col=col[6])
    T5 <- seq(0, 600, 5); lines(T5, S5(T5), col=col[7])
    T10 <- c(25, seq(100, 1000, 5)); lines(T10, S10(T10), col=col[8])
    T20 <- c(80, seq(100, 1000, 5)); lines(T20, S20(T20), col=col[9])
    T30 <- c(125, seq(200, 1000, 5)); lines(T30, S30(T30), col=col[10])
    T40 <- c(175, seq(300, 1000, 5)); lines(T40, S40(T40), col=col[11])
    T50 <- c(225, seq(300, 1000, 5)); lines(T50, S50(T50), col=col[12])
    T60 <- c(250, seq(300, 1000, 5)); lines(T60, S60(T60), col=col[13])
    legend("topleft", pch=c(0, 1, 2, 5, 6), bty = "n",
           legend=c("Helgeson, 1969", "Helgeson et al., 1981", "Manning et al., 2013", "spline control point", "high-P extrapolation"))
    legend("bottomright", col=c(NA, rev(col)), lty=1, bty = "n",
           legend=c("kbar", "60", "50", "40", "30", "20", "10", "5", "4", "3", "2", "1", "0.5", "Psat"))
    title(main=expression("Deybe-H\u00FCckel extended term ("*italic(b)[gamma]*") parameter"))
  } else if(showsplines=="P") {
    thermo.plot.new(c(0, 5), c(-.2, .7), xlab=expression(log~italic(P)*"(bar)"), ylab=expression(italic(b)[gamma]))
    # pressures that are used to make the isothermal splines (see below)
    P25 <- c(1, 500, 1000, 2000, 3000, 4000, 5000)
    P100 <- c(1, 500, 1000, 2000, 3000, 4000, 5000, 10000, 20000)
    P200 <- c(16, 500, 1000, 2000, 3000, 4000, 5000, 10000, 20000, 30000, 40000)
    P300 <- c(86, 500, 1000, 2000, 3000, 4000, 5000, 10000, 20000, 30000, 40000, 50000, 60000)
    P400 <- c(500, 1000, 2000, 3000, 4000, 5000, 10000, 20000, 30000, 40000, 50000, 60000)
    P500 <- c(1000, 2000, 3000, 4000, 5000, 10000, 20000, 30000, 40000, 50000, 60000)
    P600 <- c(2000, 3000, 4000, 5000, 10000, 20000, 30000, 40000, 50000, 60000)
    P700 <- c(10000, 20000, 30000, 40000, 50000, 60000)
    P800 <- c(10000, 20000, 30000, 40000, 50000, 60000)
    P900 <- c(10000, 20000, 30000, 40000, 50000, 60000)
    P1000 <- c(10000, 20000, 30000, 40000, 50000, 60000)
    # plot the pressure and B-dot values used to make the isothermal splines
    points(log10(P25), bgamma(25, P25))
    points(log10(P100), bgamma(100, P100))
    points(log10(P200), bgamma(200, P200))
    points(log10(P300), bgamma(300, P300))
    points(log10(P400), bgamma(400, P400))
    points(log10(P500), bgamma(500, P500))
    points(log10(P600), bgamma(600, P600))
    points(log10(P700), bgamma(700, P700))
    points(log10(P800), bgamma(800, P800))
    points(log10(P900), bgamma(900, P900))
    points(log10(P1000), bgamma(1000, P1000))
    # plot the isothermal spline functions
    col <- tail(rev(rainbow(12)), -1)
    P <- c(1, seq(50, 5000, 50)); lines(log10(P), bgamma(25, P), col=col[1])
    P <- c(1, seq(50, 20000, 50)); lines(log10(P), bgamma(100, P), col=col[2])
    P <- c(1, seq(50, 40000, 50)); lines(log10(P), bgamma(200, P), col=col[3])
    P <- c(1, seq(50, 60000, 50)); lines(log10(P), bgamma(300, P), col=col[4])
    P <- seq(500, 60000, 50); lines(log10(P), bgamma(400, P), col=col[5])
    P <- seq(1000, 60000, 50); lines(log10(P), bgamma(500, P), col=col[6])
    P <- seq(2000, 60000, 50); lines(log10(P), bgamma(600, P), col=col[7])
    P <- seq(10000, 60000, 50); lines(log10(P), bgamma(700, P), col=col[8])
    P <- seq(10000, 60000, 50); lines(log10(P), bgamma(800, P), col=col[9])
    P <- seq(10000, 60000, 50); lines(log10(P), bgamma(900, P), col=col[10])
    P <- seq(10000, 60000, 50); lines(log10(P), bgamma(1000, P), col=col[11])
    legend("topleft", col=c(NA, col), lty=1, bty = "n", legend=c("degrees C", 25, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000))
    legend("bottomright", pch=1, bty = "n", legend="points from iso-P splines")
    title(main=expression("Deybe-H\u00FCckel extended term ("*italic(b)[gamma]*") parameter"))
  } else {
    # make T and P the same length
    ncond <- max(length(T), length(P))
    T <- rep(T, length.out=ncond)
    P <- rep(P, length.out=ncond)
    # loop over P, T conditions
    bgamma <- numeric()
    lastT <- NULL
    for(i in 1:length(T)) {
      # make it fast: skip splines at 25 degC and 1 bar
      if(T[i]==25 & P[i]==1) bgamma <- c(bgamma, 0.041)
      else {
        if(!identical(T[i], lastT)) {
          # get the spline fits from particular pressures for each T
          if(T[i] >= 700) {
            PT <- c(10000, 20000, 30000, 40000, 50000, 60000)
            B <- c(S10(T[i]), S20(T[i]), S30(T[i]), S40(T[i]), S50(T[i]), S60(T[i]))
          } else if(T[i] >= 600) {
            PT <- c(2000, 3000, 4000, 5000, 10000, 20000, 30000, 40000, 50000, 60000)
            B <- c(S2(T[i]), S3(T[i]), S4(T[i]), S5(T[i]), S10(T[i]), S20(T[i]), S30(T[i]), S40(T[i]), S50(T[i]), S60(T[i]))
          } else if(T[i] >= 500) {
            PT <- c(1000, 2000, 3000, 4000, 5000, 10000, 20000, 30000, 40000, 50000, 60000)
            B <- c(S1(T[i]), S2(T[i]), S3(T[i]), S4(T[i]), S5(T[i]), S10(T[i]), S20(T[i]), S30(T[i]), S40(T[i]), S50(T[i]), S60(T[i]))
          } else if(T[i] >= 400) {
            PT <- c(500, 1000, 2000, 3000, 4000, 5000, 10000, 20000, 30000, 40000, 50000, 60000)
            B <- c(S0.5(T[i]), S1(T[i]), S2(T[i]), S3(T[i]), S4(T[i]), S5(T[i]), S10(T[i]), S20(T[i]), S30(T[i]), S40(T[i]), S50(T[i]), S60(T[i]))
          } else if(T[i] >= 300) {
            # here the lowest P is Psat
            PT <- c(86, 500, 1000, 2000, 3000, 4000, 5000, 10000, 20000, 30000, 40000, 50000, 60000)
            B <- c(S0(T[i]), S0.5(T[i]), S1(T[i]), S2(T[i]), S3(T[i]), S4(T[i]), S5(T[i]), S10(T[i]), S20(T[i]), S30(T[i]), S40(T[i]), S50(T[i]), S60(T[i]))
          } else if(T[i] >= 200) {
            # drop highest pressures because we get into ice
            PT <- c(16, 500, 1000, 2000, 3000, 4000, 5000, 10000, 20000, 30000, 40000)
            B <- c(S0(T[i]), S0.5(T[i]), S1(T[i]), S2(T[i]), S3(T[i]), S4(T[i]), S5(T[i]), S10(T[i]), S20(T[i]), S30(T[i]), S40(T[i]))
          } else if(T[i] >= 100) {
            PT <- c(1, 500, 1000, 2000, 3000, 4000, 5000, 10000, 20000)
            B <- c(S0(T[i]), S0.5(T[i]), S1(T[i]), S2(T[i]), S3(T[i]), S4(T[i]), S5(T[i]), S10(T[i]), S20(T[i]))
          } else if(T[i] >= 0) {
            PT <- c(1, 500, 1000, 2000, 3000, 4000, 5000)
            B <- c(S0(T[i]), S0.5(T[i]), S1(T[i]), S2(T[i]), S3(T[i]), S4(T[i]), S5(T[i]))
          }
          # make a new spline as a function of pressure at this T
          ST <- splinefun(PT, B)
          # remember this T; if it's the same as the next one, we won't re-make the spline
          lastT <- T[i]
        }
        bgamma <- c(bgamma, ST(P[i]))
      }
    }
    return(bgamma)
  }
}

### unexported functions ###

Bdot <- function(TC) {
  Bdot <- splinefun(c(25, 50, 100, 150, 200, 250, 300), c(0.0418, 0.0439, 0.0468, 0.0479, 0.0456, 0.0348, 0))(TC)
  Bdot[TC > 300] <- 0
  return(Bdot)
}
