# Load default settings for CHNOSZ
reset()

info <- "Simple buffer works in 0, 1, and 2 dimensions"

add.OBIGT("SUPCRT92")
# Define 4 temperatures
T <- c(200, 300, 400, 500)
# calculate SiO2 activity buffered by quartz at 1000 bar
logaSiO2 <- subcrt(c("quartz", "SiO2"), c(-1, 1), T = T, P = 1000)$out$logK

# Set up system
basis(c("Al+3", "SiO2", "Na+", "K+", "H2O", "O2", "H+"))
species(c("K-feldspar", "albite", "paragonite", "dickite", "muscovite"))
# Calculate logact(SiO2) with quartz buffer
basis("SiO2", "quartz")

# 0 dimensions: constant T
a0 <- affinity(T = 200, P = 1000)
expect_equal(a0$buffer[[1]], logaSiO2[1], info = info)

# 1 dimension: variable T
a1 <- affinity(T = T, P = 1000)
expect_equal(a1$buffer[[1]], logaSiO2, info = info)

# 2 dimensions: variable T and Na+ activity (affinity errored in version <= 1.3.6)
a2 <- affinity(T = c(200, 500, 4), "Na+" = c(1, 5, 5), P = 1000)
expect_equivalent(a2$buffer[[1]][, 1], logaSiO2, info = info)

# 2 dimensions: variable K+ and Na+ activity (diagram errored in version <= 1.3.6)
A2 <- affinity("K+" = c(1, 5, 5), "Na+" = c(1, 5, 5), T = 200, P = 1000)
#D2 <- diagram(A2)
expect_equivalent(unique(as.vector(A2$buffer[[1]])), logaSiO2[1], info = info)

# A test for 0 dimensions with a two-mineral buffer 20201103
info <- "2-mineral buffer at a single point"
# (buffer errored trying to index columns of non-matrix prior to 1.3.6-85)
T <- 400
P <- 1000
basis(c("Fe", "O2"), c("cr", "gas"))
basis("O2", "HM")
O2_HM <- affinity(T = T, P = P, return.buffer = TRUE)$O2
logfO2 <- subcrt(c("hematite", "magnetite"), c(-6, 4), T = T, P = P)$out$logK
expect_equal(O2_HM, logfO2, info = info)
