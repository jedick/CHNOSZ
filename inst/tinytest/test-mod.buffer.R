# Load default settings for CHNOSZ
reset()

# Tests added on 20260110

# Buffer + ionization: relative stabilities of E. coli sigma factors on a T-pH diagram
# (sigma factors 24, 32, 38, 54, 70, i.e. RpoE, RpoH, RpoS, RpoN, RpoD)
proteins <- c("RPOE", "RP32", "RPOS", "RP54", "RPOD")

info <- "No error for proteins as buffers"

# Define the buffer species
expect_silent(mod.buffer("sigma", paste(proteins, "ECOLI", sep = "_")), info = info)
# Get the bufferd activities of 4 basis species at 25 deg C and pH 7.4
basis("CHNOS+")
expect_silent(basis(c("CO2", "NH3", "H2S", "O2"), "sigma"), info = info)
basis("pH", 7.4)
expect_silent(logact <- affinity(return.buffer = TRUE, T = 25), info = info)

info <- "Round-trip for diagram with activities buffered by proteins"

# Set the activities of the basis species from the buffer
basis(c("CO2", "NH3", "H2S", "O2"), as.numeric(logact))
species(paste(proteins, "ECOLI", sep = "_"))
# Diagram the relative stabilities of proteins
res <- 100
a <- affinity(pH = c(5, 9, 100), T = c(20, 40, 100))
d <- diagram(a, normalize = FALSE, plot.it = FALSE)
# We should get out the T and pH given above
tp <- find.tp(d$predominant)
expect_equal(round(rev(a$vals$T)[median(tp[, 1])]), 25, info = info)
expect_equal(round(a$vals$pH[median(tp[, 2])], 1), 7.4, info = info)
