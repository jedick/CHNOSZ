# CHNOSZ/EOSregress.R  
# model volumes and heat capacities of aqueous species
# 20091105 first version
# 20110429 revise and merge with CHNOSZ package

Cp_s_var <- function(T=298.15, P=1, omega.PrTr=0, Z=0) {
  # solvation contribution to heat capacity in the HKF EOS, divided by omega(Pr,Tr) (calories)
  Cp_s <- hkf("Cp", parameters=data.frame(omega=omega.PrTr, Z=Z), T=T, P=P, contrib="s")$aq
  return(Cp_s[[1]][, 1]/omega.PrTr)
}

V_s_var <- function(T=298.15, P=1, omega.PrTr=0, Z=0) {
  # solvation contribution to volume in the HKF EOS, divided by omega(Pr,Tr) (cm3.bar)
  # [the negative sign on this term as written in the HKF EOS is accounted for by hkf()]
  V_s <- hkf("V", parameters=data.frame(omega=omega.PrTr, Z=Z), T=T, P=P, contrib="s")$aq
  return(V_s[[1]][, 1]/convert(omega.PrTr, "cm3bar"))
}

EOSvar <- function(var, T, P, ...) {
  # get the variables of a term in a regression equation
  # T (K), P (bar)
  Theta <- 228 # K
  Psi <- 2600  # bar
  out <- switch(EXPR = var,
    "(Intercept)" = rep(1, length(T)),
    "T" = T,
    "P" = P,
    "TTheta" = T - Theta,                 # T-Theta
    "invTTheta" = (T - Theta)^-1,         # 1/(T-Theta)
    "TTheta2" = (T - Theta)^2,            # (T-Theta)^2
    "invTTheta2" = (T - Theta)^-2,        # 1/(T-Theta)^2
    "invPPsi" = (P + Psi)^-1,             # 1/(P+Psi)
    "invPPsiTTheta" = (P + Psi)^-1 * (T - Theta)^-1,  # 1/[(P+Psi)(T-Theta)]
    "TXBorn" = T*water("XBorn", T=T, P=P)[, 1],
    "drho.dT" = -water("rho", T=T, P=P)[, 1]*water("E", T=T, P=P)[, 1],
    "V.kT" = water("V", T=T, P=P)[, 1]*water("kT", T=T, P=P)[, 1],
    # fallback: get a variable that is a property of water, or
    # is any other function by name (possibly a user-defined function)
    (  if(var %in% water.SUPCRT92()) water(var, T, P)[, 1]
       else if(exists(var)) {
         if(is.function(get(var))) {
           if(all(c("T", "P") %in% names(formals(get(var))))) get(var)(T=T, P=P, ...)
           else stop(paste("the arguments of ", var, "() do not contain T and P", sep=""))
         }
         else stop(paste("an object named", var, "is not a function"))
       }
       else stop(paste("can't find a variable named", var))
    )
  )
  # 20151126 apply the negative sign in the HKF EOS for V to the variable
  # (not to omega as previously assumed)
  if(var=="QBorn") out <- -out
  return(out)
}

EOSlab <- function(var, coeff="") {
  # make pretty labels for the variables
  lab <- switch(EXPR = var,
    # these are regression variables listed in EOSregress.Rd
    "(Intercept)" = substitute(YYY*" ", list(YYY=coeff)),
    "T" = substitute(YYY%*%italic(T), list(YYY=coeff)),
    "P" = substitute(YYY%*%italic(P), list(YYY=coeff)),
    "TTheta" = substitute(YYY%*%(italic(T)-Theta), list(YYY=coeff)),
    "invTTheta" = substitute(YYY/(italic(T)-Theta), list(YYY=coeff)),
    "TTheta2" = substitute(YYY%*%(italic(T)-Theta)^2, list(YYY=coeff)),
    "invTTheta2" = substitute(YYY/(italic(T)-Theta)^2, list(YYY=coeff)),
    "invPPsi" = substitute(YYY/(italic(P)+Psi),list(YYY=coeff)),
    "invPPsiTTheta" = substitute(YYY/((italic(P)+Psi)(italic(T)-Theta)), list(YYY=coeff)),
    "TXBorn" = substitute(YYY%*%italic(TX), list(YYY=coeff)),
    "drho.dT" = substitute(YYY%*%(d~rho/dT), list(YYY=coeff)),
    "V.kT" = substitute(YYY%*%V~kappa[italic(T)], list(YYY=coeff)),
    # these are non-single-letter properties of water as listed in water.Rd
    "kT" = substitute(YYY%*%kappa[italic(T)], list(YYY=coeff)),
    "alpha" = substitute(YYY%*%alpha, list(YYY=coeff)),
    "beta" = substitute(YYY%*%beta, list(YYY=coeff)),
    "epsilon" = substitute(YYY%*%epsilon, list(YYY=coeff)),
    "rho" = substitute(YYY%*%rho, list(YYY=coeff)),
    "NBorn" = substitute(YYY%*%italic(N), list(YYY=coeff)),
    "QBorn" = substitute(YYY%*%italic(Q), list(YYY=coeff)),
    "XBorn" = substitute(YYY%*%italic(X), list(YYY=coeff)),
    "YBorn" = substitute(YYY%*%italic(Y), list(YYY=coeff)),
    "ZBorn" = substitute(YYY%*%italic(Z), list(YYY=coeff)),
    (
      # if var is a function, does have an attribute named "label"?
      if(exists(var)) {
        if(is.function(get(var))) {
          if(!is.null(attr(get(var), "label"))) {
            return(substitute(YYY*XXX, list(YYY=coeff, XXX=attr(get(var), "label"))))
            # fallback, use the name of the variable
            # (e.g. for a property of water such as A, G, S, U, H, or name of a user-defined function)
          } else substitute(YYY%*%italic(XXX), list(YYY=coeff, XXX=var))
        } else substitute(YYY%*%italic(XXX), list(YYY=coeff, XXX=var))
      } else substitute(YYY%*%italic(XXX), list(YYY=coeff, XXX=var))
    )
  )
  return(lab)
}

EOSregress <- function(exptdata, var="", T.max=9999, ...) {
  # regress exptdata using terms listed in fun 
  # which values to use
  iT <- which(exptdata$T <= T.max)
  exptdata <- exptdata[iT, ]
  # temperature and pressure
  T <- exptdata$T
  P <- exptdata$P
  # the third column is the property of interest: Cp or V
  X <- exptdata[, 3]
  # now build a regression formula 
  if(length(var) == 0) stop("var is missing")
  fmla <- as.formula(paste("X ~ ", paste(var, collapse="+")))
  # retrieve the values of the variables
  for(i in seq_along(var)) assign(var[i], EOSvar(var[i], T=T, P=P, ...))
  # now regress away!
  EOSlm <- lm(fmla)
  return(EOSlm)
}

EOScalc <- function(coefficients, T, P, ...) {
  # calculate values of volume
  # or heat capacity from regression fit
  X <- 0
  for(i in 1:length(coefficients)) {
    coeff.i <- coefficients[[i]]
    fun.i <- EOSvar(names(coefficients)[i], T, P, ...)
    X <- X + coeff.i * fun.i
  }
  return(X)
}

EOSplot <- function(exptdata, var=NULL, T.max=9999, T.plot=NULL,
  fun.legend="topleft", coefficients=NULL, add=FALSE,
  lty=par("lty"), col=par("col"), ...) {
  # plot experimental and modelled volumes and heat capacities
  # first figure out the property (Cp or V) from the exptdata
  prop <- colnames(exptdata)[3]
  # if var is NULL use HKF equations
  if(is.null(var)) {
    if(prop=="Cp") var <- c("invTTheta2","TXBorn")
    if(prop=="V") var <- c("invTTheta","QBorn")
  }
  # perform the regression, only using temperatures up to T.max
  if(is.null(coefficients)) {
    EOSlm <- EOSregress(exptdata, var, T.max, ...)
    coefficients <- EOSlm$coefficients
  }
  # only plot points below a certain temperature
  iexpt <- 1:nrow(exptdata)
  if(!is.null(T.plot)) iexpt <- which(exptdata$T < T.plot)
  # for a nicer plot, extend the ranges, but don't go below -20 degrees C
  ylim <- extendrange(exptdata[iexpt, prop], f=0.1)
  xlim <- extendrange(exptdata$T[iexpt], f=0.1)
  xlim[xlim < 253.15] <- 253.15
  # start plot
  if(!add) {
    thermo.plot.new(xlim=xlim, ylim=ylim, xlab=axis.label("T", units="K"),
      ylab=axis.label(paste(prop, "0", sep="")), yline=2, mar=NULL)
    # different plot symbols to represent size of residuals
    pch.open <- 1
    pch.filled <- 16
    # find the calculated values at these conditions
    calc.X <- EOScalc(coefficients, exptdata$T, exptdata$P, ...)
    expt.X <- exptdata[, prop]
    # are we within 10% of the values
    in10 <- which(abs((calc.X-expt.X)/expt.X) < 0.1)
    pch <- rep(pch.open, length(exptdata$T))
    pch[in10] <- pch.filled
    points(exptdata$T, exptdata[, prop], pch=pch)
  }
  # plot regression line at a single P
  P <- mean(exptdata$P)
  message("EOSplot: plotting line for P=", P, " bar")
  xs <- seq(xlim[1], xlim[2], length.out=200)
  calc.X <- EOScalc(coefficients, xs, P, ...)
  lines(xs, calc.X, lty=lty, col=col)
  # make legend
  if(!is.null(fun.legend) & !add) {
    # 20161101: negate QBorn and V_s_var
    iQ <- names(coefficients) %in% c("QBorn", "V_s_var")
    coefficients[iQ] <- -coefficients[iQ]
    coeffs <- as.character(round(as.numeric(coefficients), 4))
    # so that positive ones appear with a plus sign
    ipos <- which(coeffs >= 0)
    coeffs[ipos] <- paste("+", coeffs[ipos], sep="")
    # make labels for the functions
    fun.lab <- as.expression(lapply(1:length(coeffs),
      function(x) {EOSlab(names(coefficients)[x],coeffs[x])} ))
    #fun.lab <- paste(names(coeffs),round(as.numeric(coeffs),4))
    legend(fun.legend, legend=fun.lab, pt.cex=0.1)
  }
  return(invisible(list(xrange=range(exptdata$T[iexpt]), coefficients=coefficients)))
}

EOScoeffs <- function(species, property, P=1) {
  # get the HKF coefficients for species in the database
  iis <- info(info(species, "aq"))
  if(property=="Cp") {
    out <- as.numeric(iis[,c("c1", "c2", "omega")])
    names(out) <- c("(Intercept)", "invTTheta2", "TXBorn")
  } else if(property=="V") {
    iis <- iis[,c("a1", "a2", "a3", "a4", "omega")]
    # calculate sigma and xi and convert to volumetric units: 1 cal = 41.84 cm^3 bar
    sigma <- convert( iis$a1 + iis$a2 / (2600 + P), "cm3bar" )
    xi <- convert( iis$a3 + iis$a4 / (2600 + P), "cm3bar" )
    omega <- convert( iis$omega, "cm3bar" )
    # 20151126: we _don't_ put a negative sign on omega here;
    # now, the negative sign in the HKF EOS is with the variable (QBorn or V_s_var)
    out <- c(sigma, xi, omega)
    names(out) <- c("(Intercept)", "invTTheta", "QBorn")
  }
  return(out)
}

