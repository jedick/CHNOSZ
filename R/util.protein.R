# CHNOSZ/util.protein.R
# MP90.cp: additive heat capacity from groups of Makhatadze and Privalov, 1990

MP90.cp <- function(protein, T) {
  # T (temperature, degrees C), protein (name of protein)
  # returns heat capacity of protein (kj/mol)
  # using algorithm of makhatadze and privalov, 1990.
  TMP <- c(5,25,50,75,100,125)
  A.cp <- splinefun(TMP,c(175.7,166.7,156.2,144.7,134.6,124.1))
  C.cp <- splinefun(TMP,c(225.4,237.6,250.8,260.7,268.2,276.1))
  D.cp <- splinefun(TMP,c( 72.8, 89.0,106.2,124.5,140.7,154.3))
  E.cp <- splinefun(TMP,c(168.3,179.0,192.0,203.7,211.4,217.8))
  F.cp <- splinefun(TMP,c(395.7,383.0,370.3,358.4,348.3,339.6))
  G.cp <- splinefun(TMP,c( 82.3, 78.0, 71.7, 66.4, 59.7, 53.9))
  H.cp <- splinefun(TMP,c(205.7,179.6,177.2,179.6,187.1,196.8))
  I.cp <- splinefun(TMP,c(406.8,402.3,397.1,390.8,386.0,380.8))
  K.cp <- splinefun(TMP,c(328.8,332.5,334.0,337.5,339.4,343.6))
  L.cp <- splinefun(TMP,c(385.9,381.7,377.8,372.9,369.4,365.5))
  M.cp <- splinefun(TMP,c(197.1,175.9,158.1,150.3,148.1,143.9))
  N.cp <- splinefun(TMP,c( 72.9, 88.8,109.8,125.2,140.5,154.2))
  P.cp <- splinefun(TMP,c(214.6,177.7,152.3,142.8,135.6,130.1))
  Q.cp <- splinefun(TMP,c(168.0,180.2,193.4,203.3,210.8,218.7))
  R.cp <- splinefun(TMP,c(204.6,273.4,305.8,315.1,318.7,318.5))
  S.cp <- splinefun(TMP,c( 75.6, 81.2, 85.7, 91.4, 97.3,102.1))
  T.cp <- splinefun(TMP,c(194.2,184.5,182.2,186.5,199.0,216.2))
  V.cp <- splinefun(TMP,c(324.6,314.4,305.0,294.7,285.7,269.6))
  W.cp <- splinefun(TMP,c(471.2,458.5,445.8,433.9,423.8,415.1))
  Y.cp <- splinefun(TMP,c(310.6,301.7,295.2,294.5,300.1,304.0))
  AA.cp <- splinefun(TMP,c(-158.3,-90.4,-21.5,-32.3,-92.4,-150.0))
  UPBB.cp <- splinefun(TMP,c(3.7,15.2,26.2,29.8,33.7,33.7))
  cnew <- numeric()
  for(i in 1:length(T)) {
    Ti <- T[i]
    cp <- c(A.cp(Ti),C.cp(Ti),D.cp(Ti),E.cp(Ti),F.cp(Ti),
            G.cp(Ti),H.cp(Ti),I.cp(Ti),K.cp(Ti),L.cp(Ti),
            M.cp(Ti),N.cp(Ti),P.cp(Ti),Q.cp(Ti),R.cp(Ti),
            S.cp(Ti),T.cp(Ti),V.cp(Ti),W.cp(Ti),Y.cp(Ti))
    # get the protein composition
    tt <- pinfo(pinfo(protein))[,6:25]
    cnew <- c(cnew, sum(cp * as.numeric(tt)) + sum(as.numeric(tt)) * UPBB.cp(Ti))
  }
  return(cnew)
}

### unexported functions ###

group.formulas <- function() {
  # return a matrix with chemical formulas of residues
  # names of the sidechain groups
  groups <- paste("[", aminoacids(3), "]", sep="")
  # the indices of H2O, sidechain groups, and [UPBB]
  ig <- suppressMessages(info(c("H2O", groups, "[UPBB]")))
  # their formulas
  A <- i2A(ig)
  # add [UPBB] to the sidechain groups to get residues
  out <- A[1:21,]
  out[2:21,] <- t(t(A) + A[22,])[2:21,]
  # make "H2O" not "water"
  rownames(out)[1] <- "H2O"
  return(out)
}
