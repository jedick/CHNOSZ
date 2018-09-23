# compare.R: compare SUPCRT/slop and CHNOSZ databases
# 20170228 jmd

# writes the following files:
# onot.csv - entries in thermo$obigt not in SUPCRT datafile
# snot.csv - entries in SUPCRT datafile not in thermo$obigt
# same.csv - entries that have "same" data in SUPCRT datafile and thermo$obigt (NA's not included in comparison)
# different.csv - entries that have different data in SUPCRT datafile and thermo$obigt; the numerical differences are given here

# read the SUPCRT/slop file
#sfile <- "SPRONS92.edit" # SPRONS92.DAT edited to add separator line after aq block
#sfile <- "slop98.edit"  # slop98.dat edited to add separator line, and other changes to make read.supcrt work
#sfile <- "slop07.edit"  # slop07.dat edited to remove non-ASCII characters in header (dashes in page ranges)
                       # and extra characters in ref: lines (DAPHNITE,14A and URANINITE)
sfile <- "slop15.edit"

source("read.supcrt.R")
sdat <- read.supcrt(sfile)

# put the data for cr species in the format of thermo$obigt
icr <- sdat$ghs$state=="cr"
scr <- thermo$obigt[1:sum(icr), ]
scr[] <- NA
scr$name <- sdat$ghs$name[icr]
scr$abbrv <- sdat$ghs$abbrv[icr]
scr$formula <- sdat$ghs$formula[icr]
scr$state <- sdat$ghs$state[icr]
scr$ref1 <- sdat$ghs$source[icr]
scr$date <- sdat$ghs$date[icr]
scr$G <- sdat$ghs$Gf[icr]
scr$H <- sdat$ghs$Hf[icr]
scr$S <- sdat$ghs$S[icr]
scr$a1.a <- sdat$eos.cr$a
scr$a2.b <- sdat$eos.cr$b
scr$a3.c <- sdat$eos.cr$c
scr$V <- sdat$eos.cr$V
scr$z.T <- sdat$eos.cr$T

# put the data for gas species in the format of thermo$obigt
igas <- sdat$ghs$state=="gas"
sgas <- thermo$obigt[1:sum(igas), ]
sgas[] <- NA
sgas$name <- sdat$ghs$name[igas]
sgas$abbrv <- sdat$ghs$abbrv[igas]
sgas$formula <- sdat$ghs$formula[igas]
sgas$state <- sdat$ghs$state[igas]
sgas$ref1 <- sdat$ghs$source[igas]
sgas$date <- sdat$ghs$date[igas]
sgas$G <- sdat$ghs$Gf[igas]
sgas$H <- sdat$ghs$Hf[igas]
sgas$S <- sdat$ghs$S[igas]
sgas$a1.a <- sdat$eos.gas$a
sgas$a2.b <- sdat$eos.gas$b
sgas$a3.c <- sdat$eos.gas$c
sgas$V <- sdat$eos.gas$V
sgas$z.T <- sdat$eos.gas$T

# put the data for aq species in the format of thermo$obigt
iaq <- sdat$ghs$state=="aq"
saq <- thermo$obigt[1:sum(iaq), ]
saq[] <- NA
saq$name <- sdat$ghs$name[iaq]
saq$abbrv <- sdat$ghs$abbrv[iaq]
saq$formula <- sdat$ghs$formula[iaq]
saq$state <- sdat$ghs$state[iaq]
saq$ref1 <- sdat$ghs$source[iaq]
saq$date <- sdat$ghs$date[iaq]
saq$G <- sdat$ghs$Gf[iaq]
saq$H <- sdat$ghs$Hf[iaq]
saq$S <- sdat$ghs$S[iaq]
saq$a1.a <- sdat$eos.aq$a1
saq$a2.b <- sdat$eos.aq$a2
saq$a3.c <- sdat$eos.aq$a3
saq$a4.d <- sdat$eos.aq$a4
saq$c1.e <- sdat$eos.aq$c1
saq$c2.f <- sdat$eos.aq$c2
saq$omega.lambda <- sdat$eos.aq$omega

# combine all states
sall <- rbind(scr, sgas, saq)

# replace old with new names
newnames <- read.csv("newnames.csv", as.is=TRUE)
inew <- match(sall$name, newnames$old)
sall$name[!is.na(inew)] <- newnames$new[na.omit(inew)]
# names of the species in thermo$obigt and sall
oname <- paste0(thermo$obigt$name, "_", thermo$obigt$state)
sname <- paste0(sall$name, "_", sall$state)
# index of matching and non-matching species in thermo$obigt and sall
omatch <- onot <- numeric()
smatch <- snot <- numeric()
# loop over species in sall
for(i in 1:nrow(sall)) {
  # try to find a match in thermo$obigt
  imatch <- match(sname[i], oname)
  # if it's NA, try cr1
  if(is.na(imatch)) imatch <- match(paste0(sname[i], "1"), oname)
  # update the indices
  if(is.na(imatch)) snot <- c(snot, i)
  else {
    smatch <- c(smatch, i)
    omatch <- c(omatch, imatch)
  }
}
# that leaves onot
onot <- (1:nrow(thermo$obigt))[-omatch]

# write the non-matching entries first
write.csv(thermo$obigt[onot, ], "onot.csv", quote=c(1, 2))
write.csv(sall[snot, ], "snot.csv", quote=c(1, 2, 5))

# calculate the difference between matching entries
oout <- thermo$obigt[omatch, ]
sout <- sall[smatch, ]
# combine ref1 and ref2
oout$ref2[is.na(oout$ref2)] <- ""
oout$ref2 <- paste0(oout$ref1, ",", oout$ref2)
oout$ref1 <- sout$ref1
# the columns with numeric data
icol <- 8:20
for(i in icol) {
  # put the difference in oout
  oout[, i] <- round(oout[, i] - as.numeric(sout[, i]), 4)
}
# find whether each entry has the same, almost the same, or different data
ialmost <- isame <- numeric()
for(i in 1:nrow(oout)) {
  issame <- all(sapply(sapply(na.omit(as.numeric(oout[i, icol])), all.equal, 0), isTRUE))
  if(issame) isame <- c(isame, i)
  # G, H, S, Cp, V, a1, a2, a3, a4, c1, c2, omega, T
  isalmost <- all(sapply(na.omit(abs(as.numeric(oout[i, icol])) < c(500, 500, 2, 1, 1, 0.1, 1, 1, 1, 0.5, 0.5, .1, 1)), isTRUE))
  if(isalmost) ialmost <- c(ialmost, i)
}
# write the difference between matching entries
write.csv(oout[isame, ], "same.csv", quote=c(1, 2, 5, 6))
write.csv(oout[setdiff(ialmost, isame), ], "almost.csv", quote=c(1, 2, 5, 6))
write.csv(oout[-ialmost, ], "diffs.csv", quote=c(1, 2, 5, 6))
# write the values from the SUPCRT file for the different entries
write.csv(sout[-ialmost, ], "different.csv", quote=c(1, 2, 5, 6))
