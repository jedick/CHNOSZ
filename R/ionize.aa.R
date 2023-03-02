# CHNOSZ/ionize.aa.R
# Rewritten ionization function 20120526 jmd

ionize.aa <- function(aa, property="Z", T=25, P="Psat", pH=7, ret.val=NULL, suppress.Cys=FALSE) {
  # Calculate the additive ionization property of proteins with amino acid 
  # composition in aa as a function of vectors of T, P and pH;
  # property if NULL is the net charge, if not NULL is one of the subcrt() properties
  # or "A" to calculate A/2.303RT for the ionization reaction
  # ret.val can be 'pK', 'alpha' or 'aavals' to get these values;
  # T, P and pH should be same length, or larger a multiple of the smaller
  lmax <- max(c(length(T), length(P), length(pH)))
  T <- rep(T, length.out=lmax)
  P <- rep(P, length.out=lmax)
  pH <- rep(pH, length.out=lmax)
  # Turn pH into a matrix with as many columns as ionizable groups
  pH <- matrix(rep(pH, 9), ncol=9)
  # Turn charges into a matrix with as many rows as T,P,pH conditions
  charges <- c(-1, -1, -1, 1, 1, 1, -1, 1, -1)
  charges <- matrix(rep(charges, lmax), nrow=lmax, byrow=TRUE)
  # The rownumbers of the ionizable groups in thermo()$OBIGT
  neutral <- c("[Cys]", "[Asp]", "[Glu]", "[His]", "[Lys]", "[Arg]", "[Tyr]", "[AABB]", "[AABB]")
  charged <- c("[Cys-]","[Asp-]","[Glu-]","[His+]","[Lys+]","[Arg+]","[Tyr-]","[AABB+]","[AABB-]")
  ineutral <- info(neutral, "aq")
  icharged <- info(charged, "aq")
  # We'll only call subcrt() with the unique pressure/temperature combinations
  pTP <- paste(T, P)
  dupPT <- duplicated(pTP)
  # What property are we after
  sprop <- c("G", property)
  if(property %in%  c("A", "Z")) sprop <- "G"
  # Use convert=FALSE so we get results in Joules 20210407
  # (means we need to supply temperature in Kelvin)
  TK <- convert(T, "K")
  sout <- subcrt(c(ineutral, icharged), T=TK[!dupPT], P=P[!dupPT], property=sprop, convert = FALSE)$out
  # The G-values
  Gs <- sapply(sout, function(x) x$G)
  # Keep it as a matrix even if we have only one unique T, P-combo
  if(length(pTP[!dupPT])==1) Gs <- t(Gs)
  # Now the Gibbs energy difference for each group
  DG <- Gs[, 10:18, drop=FALSE] - Gs[, 1:9, drop=FALSE]
  # Build a matrix with one row for each of the (possibly duplicated) T, P values
  uPT <- unique(pTP)
  DG <- t(sapply(pTP, function(x) DG[match(x, uPT), , drop=FALSE]))
  # The pK values (-logK) 
  DG <- DG * charges
  pK <- apply(DG, 2, function(x) convert(x, "logK", T = TK))
  # Keep it as a matrix even if we have only one T, P-combo
  if(lmax==1) pK <- t(pK)
  if(identical(ret.val, "pK")) {
    colnames(pK) <- charged
    return(pK)
  }
  # Now to calculate alpha! - degrees of formation of the charged groups
  alpha <- 1 / (1 + 10 ^ (charges * (pH - pK)))
  # Suppress cysteine ionization if requested
  if(suppress.Cys) alpha[, 1] <- 0
  if(identical(ret.val, "alpha")) return(alpha)
  # Now to calculate the properties of the ionizable groups - can be charges, 
  # the chemical affinities of the ionization reactions,
  # or another property from subcrt()
  if(identical(property, "Z")) aavals <- charges
  else if(identical(property, "A")) aavals <- - charges * (pH - pK)
  else {
    # It's not charge, so compile it from the subcrt output
    # the property-values
    icol <- match(property, colnames(sout[[1]]))
    aavals <- sapply(sout, function(x) x[,icol])
    # Keep it as a matrix even if we have only one unique T, P-combo
    if(length(pTP[!dupPT])==1) aavals <- t(aavals)
    # Build a matrix with one row for each of the (possibly duplicated) T, P values
    aavals <- t(sapply(pTP, function(x) aavals[match(x, uPT), , drop=FALSE], USE.NAMES=FALSE))
    # The property difference for each group
    aavals <- aavals[, 10:18, drop=FALSE] - aavals[, 1:9, drop=FALSE]
  }
  if(identical(ret.val, "aavals")) {
    colnames(aavals) <- charged
    return(aavals)
  }
  # The contribution from each group to the ionization property of the protein
  aavals <- aavals * alpha
  # Now work with 'aa'; the next line is so that a missing argument shows
  # The name of this function in the error message
  aa <- aa
  # The columns where we find the counts of ionizable groups
  iionize <- match(c("Cys", "Asp", "Glu", "His", "Lys", "Arg", "Tyr", "chains", "chains"), colnames(aa))
  aa <- as.matrix(aa[, iionize])
  # Add it all up
  out <- apply(aa, 1, function(x) {
    aavals %*% x
  })
  # Keep it as a matrix even if we have only one T, P-combo
  if(lmax==1) out <- t(out)
  rownames(out) <- rownames(aavals)
  # That's all folks!
  return(out)
}
