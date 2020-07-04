# CHNOSZ/water.R
# calculate thermodynamic and electrostatic properties of H2O
# 20061016 jmd

water <- function(property = NULL, T = 298.15, P = "Psat", P1 = TRUE) {
  # calculate the properties of liquid H2O as a function of T and P
  # T in Kelvin, P in bar
  if(is.null(property)) return(get("thermo", CHNOSZ)$opt$water)
  # set water option
  if(length(property)==1 & any(property %in% c("SUPCRT", "SUPCRT92", "IAPWS", "IAPWS95", "DEW"))) {
    # change references 20200629
    if(property %in% c("SUPCRT", "SUPCRT92")) {
      suppressMessages(mod.OBIGT("water", ref1 = "HGK84"))
      suppressMessages(mod.OBIGT("water", ref2 = "JOH92"))
    } else if(property %in% c("IAPWS", "IAPWS95")) {
      suppressMessages(mod.OBIGT("water", ref1 = "WP02"))
      suppressMessages(mod.OBIGT("water", ref2 = NA))
    } else if(property == "DEW") {
      suppressMessages(mod.OBIGT("water", ref1 = "SHA14"))
      suppressMessages(mod.OBIGT("water", ref2 = NA))
    }
    oldwat <- thermo()$opt$water
    thermo("opt$water" = property)
    message(paste("water: setting water model to", property))
    return(invisible(oldwat))
  }
  # make T and P equal length
  if(!identical(P, "Psat")) {
    if(length(P) < length(T)) P <- rep(P, length.out = length(T))
    else if(length(T) < length(P)) T <- rep(T, length.out = length(P))
  }
  wopt <- get("thermo", CHNOSZ)$opt$water
  if(grepl("SUPCRT", wopt)) {
    # change 273.15 K to 273.16 K (needed for water.SUPCRT92 at Psat)
    if(identical(P, "Psat")) T[T == 273.15] <- 273.16
    # get properties using SUPCRT92
    w.out <- water.SUPCRT92(property, T, P, P1)
  }
  if(grepl("IAPWS", wopt)) {
    # get properties using IAPWS-95 
    w.out <- water.IAPWS95(property, T, P)
  }
  if(grepl("DEW", wopt)) {
    # use the Deep Earth Water (DEW) model
    w.out <- water.DEW(property, T, P)
  }
  w.out
}

water.SUPCRT92 <- function(property=NULL, T=298.15, P=1, P1=TRUE) {
  ### interface to H2O92D.f : FORTRAN subroutine taken from 
  ### SUPCRT92 for calculating the thermodynamic and 
  ### electrostatic properties of H2O. 
  ## we restrict the calculations to liquid water
  ## except for getting Psat (vapor-liquid saturation 
  ## pressure as a function of T>100 C). 20071213 jmd
  available_properties <- c("A", "G", "S", "U", "H", "Cv", "Cp",
    "Speed", "alpha", "beta", "epsilon", "visc",
    "tcond", "surten", "tdiff", "Prndtl", "visck", "albe",
    "ZBorn", "YBorn", "QBorn", "daldT", "XBorn",
    "V", "rho", "Psat", "E", "kT",
    "A_DH", "B_DH")
  if(is.null(property)) return(available_properties)
  # check for availability of properties
  iprop <- match(property, available_properties)
  if(any(is.na(iprop))) stop(paste("property(s) not available:", paste(property[is.na(iprop)], collapse=" ")))
  # make sure Psat in properties comes with isat=1
  if("Psat" %in% property & !identical(P, "Psat")) stop("please set P='Psat' to calculate the property Psat")
  # for Psat(T) (1) or T-P (2)
  if(identical(P, "Psat")) iopt <- 1 else iopt <- 2
  if(identical(P, "Psat")) isat <- 1 else isat <- 0
  # input values, gleaned from H2O92D.f and SUP92D.f
  # it, id, ip, ih, itripl, isat, iopt, useLVS, epseqn, icrit
  specs <- c(2, 2, 2, 5, 1, isat, iopt, 1, 4, 0)
  states <- rep(0, 4)
  # initialize the output matrix
  w.out <- matrix(NA, nrow=length(T), ncol=23, byrow=TRUE) 
  err.out <- numeric(length(T))
  rho.out <- numeric(length(T))
  P.out <- numeric(length(T))
  # 20091022 TODO: parallelize this
  Tc <- convert(T, "C")
  for(i in 1:length(T)) {
    states[1] <- Tc[i]
    if(identical(P, "Psat")) states[2] <- 0
    else states[2] <- P[i]
    if(is.na(Tc[i]) | is.na(P[i]) & !identical(P, "Psat")) {
      # if T or P is NA, all properties are NA
      # (NA's are already in w.out)
      P.out[i] <- NA
      rho.out[i] <- NA
    } else {
      # now to the actual calculations
      H2O <- .Fortran(C_h2o92, as.integer(specs), as.double(states),
        as.double(rep(0, 46)), as.integer(0))
      # errors
      err.out[i] <- H2O[[4]]
      # density of two states
      rho <- H2O[[2]][3]
      rho2 <- H2O[[2]][4]
      if(rho2 > rho) {
        # liquid is denser than vapor
        rho <- rho2 
        inc <- 1  # second state is liquid
      } else inc <- 0  # first state is liquid
      rho.out[i] <- rho
      # 23 properties of the phase in the liquid state
      w <- t(H2O[[3]][seq(1, 45, length.out=23)+inc])
      if(err.out[i]==1) w[1, ] <- NA
      # update the ith row of the output matrix
      w.out[i,] <- w
      # Psat
      if(identical(P, "Psat")) {
        w.P <- H2O[[2]][2]
        w.P[w.P==0] <- NA
        # Psat specifies P=1 below 100 degC
        if(P1) w.P[w.P < 1] <- 1
        P.out[i] <- w.P
      }
    }
  }
  # convert output to dataframe
  w.out <- as.data.frame(w.out)
  # add names of properties to the output
  names(w.out) <- available_properties[1:23]
  # assemble additional properties: V, rho, Psat, E, kT
  if(any(iprop > 23)) {
    mwH2O <- 18.0152 # SUP92.f
    V <- mwH2O/rho.out
    rho <- rho.out*1000
    # rho==0 should be NA 20180923
    rho[rho==0] <- NA
    Psat <- P.out
    E <- V*w.out$alpha
    kT <- V*w.out$beta
    # A and B parameters in Debye-Huckel equation:
    # Helgeson (1969) doi:10.2475/ajs.267.7.729
    # Manning (2013) doi:10.2138/rmg.2013.76.5
    A_DH <- 1.8246e6 * rho.out^0.5 / (w.out$epsilon * T)^1.5
    B_DH <- 50.29e8 * rho.out^0.5 / (w.out$epsilon * T)^0.5
    w.out <- cbind(w.out, data.frame(V=V, rho=rho, Psat=Psat, E=E, kT=kT, A_DH=A_DH, B_DH=B_DH))
  }
  # tell the user about any problems
  if(any(err.out==1)) {
    if(length(T) > 1) plural <- "s" else plural <- ""
    nerr <- sum(err.out==1)
    if(nerr > 1) plural2 <- "s" else plural2 <- ""
    if(identical(P, "Psat")) message(paste("water.SUPCRT92: error", plural2, " calculating ",
      nerr, " of ", length(T), " point", plural, "; for Psat we need 273.16 < T < 647.067 K", sep=""))
    else message(paste("water.SUPCRT92: error", plural2, " calculating ", nerr,
      " of ", length(T), " point", plural,
      "; T < Tfusion@P, T > 2250 degC, or P > 30kb.", sep=""))
      # that last bit is taken from SUP92D.f in SUPCRT92
  }
  # return only the selected properties
  return(w.out[, iprop, drop=FALSE])
}

water.IAPWS95 <- function(property=NULL, T=298.15, P=1) {
  available_properties <- c("A", "G", "S", "U", "H", "Cv", "Cp",
    "Speed", "epsilon",
    "YBorn", "QBorn", "XBorn", "NBorn", "UBorn",
    "V", "rho", "Psat", "de.dT", "de.dP", "pressure",
    "A_DH", "B_DH")
  if(is.null(property)) return(available_properties) 
  # to get the properties of water via IAPWS-95
  message(paste("water.IAPWS95: calculating", length(T), "values for"), appendLF=FALSE)
  M <- 18.015268 # g mol-1
  V <- function() return(M*1000/my.rho)
  # Psat stuff
  Psat <- function() {
    P <- WP02.auxiliary("P.sigma", T)
    P[T < 373.124] <- 0.1
    return(convert(P, "bar"))
  }
  ## thermodynamic properties
  Tr <- 298.15
  # convert to SUPCRT reference state
  # at the triple point
  # I2S = SUPCRT - IAPWS ( + entropy in G )
  dH <- -68316.76 - 451.75437
  dS <- 16.7123 - 1.581072
  dG <- -56687.71 + 19.64228 - dS * (T - Tr)
  # does the reference state used for GHS also go here?
  dU <- -67434.5 - 451.3229
  dA <- -55814.06 + 20.07376 - dS * (T - Tr)
  # calculate pressure from the given T and estimated rho
  pressure <- function() return(convert(IAPWS95("p", T=T, rho=my.rho), "bar"))
  # convert IAPWS95() (specific, joule) to (molar, cal) 
  S <- function()
    return(convert(IAPWS95('s',T=T,rho=my.rho)$s*M,'cal')+dS) 
  # u (internal energy) is not here because the letter
  # is used to denote one of the Born functions
  # scratch that! let's put u here and call the other one uborn
  U <- function()
    return(convert(IAPWS95('u',T=T,rho=my.rho)$u*M,'cal')+dU)
  A <- function()
    return(convert(IAPWS95('a',T=T,rho=my.rho)$a*M,'cal')+dA)
  H <- function() 
    return(convert(IAPWS95('h',T=T,rho=my.rho)$h*M,'cal')+dH) 
  G <- function() 
    return(convert(IAPWS95('g',T=T,rho=my.rho)$g*M,'cal')+dG) 
  Cv <- function() 
    return(convert(IAPWS95('cv',T=T,rho=my.rho)$cv*M,'cal')) 
  Cp <- function() 
    return(convert(IAPWS95('cp',T=T,rho=my.rho)$cp*M,'cal')) 
  Speed <- function()
    return(IAPWS95('w',T=T,rho=my.rho)$w*100) # to cm/s
  ## electrostatic properties
  epsilon <- function() return(water.AW90(T=T,rho=my.rho,P=convert(P,'MPa')))
  de.dT <- function() {
    p <- numeric()
    for(i in 1:length(T)) {
      this.T <- T[i]
      this.P <- P[i]
      this.rho <- my.rho[i]
      dt <- 0.001; t1 <- this.T-dt; t2 <- this.T+dt
      rho <- rho.IAPWS95(T=c(t1, t2), P=rep(this.P, 2))
      e <- water.AW90(T=c(t1,t2),rho=rho,rep(this.P,2))
      p <- c(p,(e[2]-e[1])/(2*dt))
    }
    return(p)
  }
  de.dP <- function() {
    p <- numeric()
    for(i in 1:length(T)) {
      this.T <- T[i]
      this.P <- P[i]
      this.rho <- my.rho[i]
      dp <- 0.001; p1 <- this.P-dp; p2 <- this.P+dp
      rho <- rho.IAPWS95(T=rep(this.T, 2), P=c(p1, p2))
      e <- water.AW90(P=c(p1,p2),rho=rho,T=rep(this.T,2))
      p <- c(p,(e[2]-e[1])/(2*dp))
    }
    return(p)
  }
  ## Born functions
  QBorn <- function() {
    p <- numeric()
    for(i in 1:length(T)) {
      this.T <- T[i]; this.P <- P[i]; this.rho <- my.rho[i]
      dp <- 0.01; p1 <- this.P-dp; p2 <- this.P+dp
      rho <- rho.IAPWS95(T=rep(this.T, 2), P=c(p1, p2))
      e <- water.AW90(T=rep(this.T,2),rho=rho,P=convert(c(p1,p2),'MPa'))
      #p <- c(p,convert(-(1/e[2]-1/e[1])/(2*dp),'cm3bar'))
      p <- c(p,-(1/e[2]-1/e[1])/(2*dp))
    }
    return(p)
  }
  NBorn <- function() {
    p <- numeric()
    for(i in 1:length(T)) {
      this.T <- T[i]; this.P <- P[i]; this.rho <- my.rho[i]
      dp <- 0.01; p1 <- this.P-dp; p2 <- this.P+dp
      rho <- rho.IAPWS95(T=rep(this.T, 3), P=c(p1, this.P, p2))
      e <- water.AW90(T=rep(this.T,3),rho=rho,P=convert(c(p1,this.P,p2),'MPa'))
      #p <- c(p,convert(convert((-(1/e[3]-1/e[2])/dp+(1/e[2]-1/e[1])/dp)/dp,'cm3bar'),'cm3bar'))
      p <- c(p,(-(1/e[3]-1/e[2])/dp+(1/e[2]-1/e[1])/dp)/dp)
    }
    return(p)
  }
  YBorn <- function() {
    p <- numeric()
    for(i in 1:length(T)) {
      this.T <- T[i]; this.P <- P[i]; this.rho <- my.rho[i]
      dt <- 0.001; t1 <- this.T-dt; t2 <- this.T+dt
      rho <- rho.IAPWS95(T=c(t1, t2), P=rep(this.P, 2))
      e <- water.AW90(T=c(t1,t2),rho=rho,P=convert(rep(this.P,2),'MPa'))
      p <- c(p,-(1/e[2]-1/e[1])/(2*dt))
    }
    return(p)
  }
  XBorn <- function() {
    p <- numeric()
    for(i in 1:length(T)) {
      this.T <- T[i]; this.P <- P[i]; this.rho <- my.rho[i]
      dt <- 0.001; t1 <- this.T-dt; t2 <- this.T+dt
      rho <- rho.IAPWS95(T=c(t1, this.T, t2), P=rep(this.P, 3))
      e <- water.AW90(T=c(t1,this.T,t2),rho=rho,P=convert(rep(this.P,3),'MPa'))
      p <- c(p,(-(1/e[3]-1/e[2])/dt+(1/e[2]-1/e[1])/dt)/dt)
    }
    return(p)
  }
  UBorn <- function() {
    p <- numeric()
    for(i in 1:length(T)) {
      this.T <- T[i]; this.P <- P[i]; this.rho <- my.rho[i]
      dt <- 0.001; this.T1 <- this.T - dt; this.T2 <- this.T + dt
      dp <- 0.001; p1 <- this.P-dp; p2 <- this.P+dp
      rho1 <- rho.IAPWS95(T=rep(this.T1, 2), P=c(p1, p2))
      rho2 <- rho.IAPWS95(T=rep(this.T2, 2), P=c(p1, p2))
      e1 <- water.AW90(T=rep(this.T1,2),rho=rho1,P=convert(c(p1,p2),'MPa'))
      e2 <- water.AW90(T=rep(this.T2,2),rho=rho2,P=convert(c(p1,p2),'MPa'))
      #p1 <- convert(-(1/e1[2]-1/e1[1])/(2*dp),'cm3bar')
      #p2 <- convert(-(1/e2[2]-1/e2[1])/(2*dp),'cm3bar')
      p1 <- -(1/e1[2]-1/e1[1])/(2*dp)
      p2 <- -(1/e2[2]-1/e2[1])/(2*dp)
      p <- c(p,(p2-p1)/(2*dt))
    }
    return(p)
  }
  ### main loop; init dataframe output and density holders
  w.out <- data.frame(matrix(nrow=length(T), ncol=length(property)))
  my.rho <- NULL
  # get densities unless only Psat is requested
  if(!identical(property, "Psat")) {
    # calculate values of P for Psat
    if(identical(P, "Psat")) P <- Psat()
    message(" rho", appendLF=FALSE)
    my.rho <- rho.IAPWS95(T, P, get("thermo", CHNOSZ)$opt$IAPWS.sat)
    rho <- function() my.rho
  }
  # get epsilon and A_DH, B_DH (so we calculate epsilon only once)
  if(any(property %in% c("epsilon", "A_DH", "B_DH"))) {
    my.epsilon <- epsilon()
    epsilon <- function() my.epsilon
    A_DH <- function() 1.8246e6 * (my.rho/1000)^0.5 / (my.epsilon * T)^1.5
    B_DH <- function() 50.29e8 * (my.rho/1000)^0.5 / (my.epsilon * T)^0.5
  }
  for(i in 1:length(property)) {
    if(property[i] %in% c("E", "kT", "alpha", "daldT", "beta")) {
      # E and kT aren't here yet... set them to NA
      # also set alpha, daldT, and beta (for derivatives of g function) to NA 20170926 
      warning("water.IAPWS95: values of ", property[i], " are NA", call.=FALSE)
      wnew <- rep(NA, length(T))
    } else {
      message(paste(" ", property[i], sep=""), appendLF=FALSE)
      wnew <- get(property[i])()
    }
    w.out[, i] <- wnew
  }  
  message("")
  # include properties available in SUPCRT that might be NA here
  wprop <- unique(c(water.SUPCRT92(), available_properties))
  iprop <- match(property, wprop)
  property[!is.na(iprop)] <- wprop[na.omit(iprop)]
  colnames(w.out) <- property
  return(w.out)
}

# get water properties from DEW model for use by subcrt() 20170925
water.DEW <- function(property = NULL, T = 373.15, P = 1000) {
  available_properties <- c("G", "epsilon", "QBorn", "V", "rho", "beta", "A_DH", "B_DH")
  if(is.null(property)) return(available_properties)
  # we can't use Psat here
  if(identical(P, "Psat")) stop("Psat isn't supported in this implementation of the DEW model. Try setting P to at least 1000 bar.")
  # use uppercase property names (including H, S, etc., so we get them from the SUPCRT properties)
  wprop <- water.SUPCRT92()
  iprop <- match(property, wprop)
  property[!is.na(iprop)] <- wprop[na.omit(iprop)]
  # convert temperature units
  pressure <- P
  temperature <- convert(T, "C")
  # initialize output data frame with NA for all properties and conditions
  ncond <- max(length(T), length(P))
  out <- matrix(NA, ncol=length(property), nrow=ncond)
  out <- as.data.frame(out)
  colnames(out) <- property
  # calculate rho and epsilon if they're needed for any other properties
  if(any(c("rho", "V", "QBorn", "epsilon", "beta", "A_DH", "B_DH") %in% property)) rho <- calculateDensity(pressure, temperature)
  if(any(c("epsilon", "A_DH", "B_DH") %in% property)) epsilon <- calculateEpsilon(rho, temperature)
  # fill in columns with values
  if("rho" %in% property) out$rho <- rho*1000 # use kg/m^3 (like SUPCRT)
  if("V" %in% property) out$V <- 18.01528/rho
  if("G" %in% property) out$G <- calculateGibbsOfWater(pressure, temperature)
  if("QBorn" %in% property) out$QBorn <- calculateQ(rho, temperature)
  if("epsilon" %in% property) out$epsilon <- epsilon
  # divide drhodP by rho to get units of bar^-1
  if("beta" %in% property) out$beta <- calculate_drhodP(rho, temperature) / rho
  if("A_DH" %in% property) out$A_DH <- 1.8246e6 * rho^0.5 / (epsilon * T)^1.5
  if("B_DH" %in% property) out$B_DH <- 50.29e8 * rho^0.5 / (epsilon * T)^0.5
  # use SUPCRT-calculated values below 100 degC and/or below 1000 bar
  ilow <- T < 373.15 | P < 1000
  if(any(ilow)) {
    out[ilow, ] <- water.SUPCRT92(property, T=T[ilow], P=P[ilow])
    iPrTr <- T == 298.15 & P == 1
    if(sum(iPrTr)==sum(ilow)) message(paste("water.DEW: using SUPCRT calculations for Pr,Tr"))
    if(sum(iPrTr)==0) message(paste("water.DEW: using SUPCRT calculations for", sum(ilow), "low-T or low-P condition(s)"))
    if(sum(iPrTr)==1 & sum(ilow) > sum(iPrTr)) message(paste("water.DEW: using SUPCRT calculations for Pr,Tr and", sum(ilow)-1, "other low-T or low-P condition(s)"))
    # however, we also have this:
    # epsilon(Pr,Tr) from SUPCRT: 78.24514
    # epsilon(Pr,Tr) in DEW spreadsheet: 78.47
    if("epsilon" %in% property) out$epsilon[iPrTr] <- 78.47
  }
  out
}
