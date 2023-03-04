# CHNOSZ/demo/carboxylase.R
# Animate rank-activity diagrams over a temperature
# and logaH2 gradient, or plot a single one for a single temperature
# First version ca. 200903; packaged anim.carboxylase() 20110818; converted to demo 20171030
library(CHNOSZ)

T <- 25:125
ntop <- 5
lcex <- 0.8
width <- 380
height <- 380
# Plot rank-activity diagrams for 24 carboxylases;
# 12 ribulose phosphate carboxylase
# 12 acetyl-coenzyme A carboxylase
# 6 of each type are nominally from mesophilic organisms
# and 6 from thermophilic organisms
# arranged here in order of increasing growth temperature
rubisco <- c("RBL_BRAJA", "A6YF84_9PROT", "A1E8R4_9CHLO", "A8C9T6_9MYCO", "A3EQE1_9BACT", "A5CKC7_9CHRO", 
  "RBL_SYNJA", "Q6JAI0_9RHOD", "RBL_METJA", "A3DND9_STAMF", "A1RZJ5_THEPD", "RBL_PYRHO")
rubisco.organisms <- c("a-proteobacterium-R", "b-proteobacterium", "Bracteacoccus", "Mycobacterium", 
  "Leptospirillum", "Cyanobium", "Synechococcus", "Cyanidiales", 
  "Methanococcus-R", "Desulfurococcus", "Thermofilum", "Pyrococcus")
accoaco <- c("Q9F7M8_PRB01", "ACCA_DEIRA", "A6CDM2_9PLAN", "A4AGS7_9ACTN", "ACCA_CAUCR", "A1VC70_DESVV", 
  "A6VIX9_METM7", "Q2JSS7_SYNJA", "A0GZU2_9CHLR", "A7WGI1_9AQUI", "Q05KD0_HYDTH", "ACCA_AQUAE")
accoaco.organisms <- c("g-proteobacterium", "Deinococcus", "Planctomyces", "Actinobacterium", 
  "a-proteobacterium-A", "d-proteobacterium", "Methanococcus-A", "Synechococcus", 
  "Chloroflexus", "Hydrogenobaculum", "Hydrogenobacter", "Aquifex")
# Assemble them all
organisms <- c(rubisco.organisms, accoaco.organisms)
# New scheme 20090611: red for hot,  blue for cold
# Open for rubisco,  filled for accoaco
col <- rep(c(rep("blue", 6), rep("red", 6)), 2)
pch <- c(rep(c(0:2, 5:7), 2), rep(c(15:20), 2))
# How many frames do we want?
res <- length(T)
if(res == 1) ido <- 1 else {
  # Check for png directory
  if(!"png" %in% dir()) stop("directory 'png' not present")
  else if(length(dir("png")) > 0) stop("directory 'png' not empty")
  # Start the plot device - multiple png figures
  png(filename = "png/Rplot%04d.png", width = width, height = height)
  # Add counters for lead-in and lead-out frames
  ido <- c(rep(1, 15), 1:res, rep(res, 20))
}
# Set up system
basis(c("CO2", "H2O", "NH3", "H2", "H2S", "H+"), 
  c("aq", "liq", "aq", "aq", "aq", "aq"), c(-3, 0, -4, -6, -7, -7))
species(c(rubisco,accoaco))
# Equation for logaH2 as a function of temperature
# from Dick and Shock, 2011
# http://dx.plos.org/10.1371/journal.pone.0022782
get.logaH2 <- function(T) return(-11 + T * 3 / 40)
H2 <- get.logaH2(T)
# Calculate affinities
if(res == 1) {
  basis("H2",H2)
  a <- affinity(T = T)
} else a <- affinity(T = T,H2 = H2)
# Calculate activities
e <- equilibrate(a, normalize = TRUE)
# For each point make a rank plot
rank <- 1:length(e$loga.equil)
for(i in 1:length(ido)) {
  # Print some progress
  if(i%%20 == 0) cat("\n") else cat(".")
  # Keep track of positions of previous points
  loga <- numeric()
  for(j in 1:length(e$loga.equil)) loga <- c(loga, e$loga.equil[[j]][ido[i]])
  if(i > 4) myrank4 <- myrank3
  if(i > 3) myrank3 <- myrank2
  if(i > 2) myrank2 <- myrank1
  if(i > 1) myrank1 <- myrank
  order <- order(loga,decreasing = TRUE)
  myrank <- rank(loga)
  cex <- rep(1.2,24)
  # Show changes by increasing point size
  # Any points that changed on the step before the step before the step before?
  if(i > 4) {
    ichanged <- myrank3 != myrank4
    cex[ichanged[order]] <- cex[ichanged[order]] + 0.1
  }
  # Any points that changed on the step before the step before?
  if(i > 3) {
    ichanged <- myrank2 != myrank3
    cex[ichanged[order]] <- cex[ichanged[order]] + 0.2
  }
  # Any points that changed on the step before?
  if(i > 2) {
    ichanged <- myrank1 != myrank2
    cex[ichanged[order]] <- cex[ichanged[order]] + 0.3
  }
  # Any points that changed on this step?
  if(i > 1) {
    ichanged <- myrank != myrank1
    cex[ichanged[order]] <- cex[ichanged[order]] + 0.4
  }
  plot(rank,loga[order],col = col[order],pch = pch[order],
    ylab = expression(log~italic(a)), cex = cex, cex.main = 1, cex.lab = 1, cex.axis = 1)
  myT <- format(round(T, 1))[ido[i]]
  myH2 <- format(round(H2, 2))[ido[i]]
  title(main = substitute(list(X~degree*C, log*italic(a)[paste(H2)] == Y), 
    list(X = myT, Y = myH2)))
  # Legends showing highest and lowest few
  legend("topright", legend = c(paste("top", ntop), organisms[order[1:ntop]]), 
    pch = c(NA, pch[order[1:ntop]]), col = c(NA, col[order[1:ntop]]), 
    pt.cex = c(NA, cex[1:ntop]), cex = lcex)
  order <- order(loga)
  legend("bottomleft", legend = c(paste("low", ntop), organisms[order[ntop:1]]), 
    pch = c(NA, pch[order[ntop:1]]), col = c(NA, col[order[ntop:1]]), 
    pt.cex = c(NA, cex[24:(24-ntop+1)]), cex = lcex)
}
# Finish up animation stuff
if(res > 1) {
  # Finish progress report
  cat("\n")
  # Close PNG plot device
  dev.off()
  # Make animated GIF using ImageMagick
  cat("anim.carboxylase: converting to animated GIF...\n")
  outfile <- "carboxylase.gif"
  syscmd <- paste("convert -loop 0 -delay 10 png/*.png png/", outfile, sep = "")
  cat(paste(syscmd,"\n"))
  if(.Platform$OS.type == "unix") sres <- system(syscmd)
  else sres <- shell(syscmd)
  if(sres == 0) cat(paste("anim.carboxylase: animation is at png/", outfile, "\n", sep = ""))
  else {
    cat("anim.carboxylase: error converting to animated GIF\n")
    cat("anim.carboxylase: check that 'convert' tool from ImageMagick is in your PATH\n")
  }
}
