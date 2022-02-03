# test-berman.R 20171001

# Load default settings for CHNOSZ
reset()

# Make sure all Berman minerals are listed with units of J in OBIGT 20220203
info <- "Berman minerals are listed with units of J in OBIGT"
file <- system.file("extdata/OBIGT/Berman_cr.csv", package = "CHNOSZ")
dat <- read.csv(file)
expect_true(all(dat$E_units == "J"), info = info)

# The maximum absolute pairwise difference between x and y
maxdiff <- function(x, y) max(abs(y - x))

info <- "high-T,P calculated properties are similar to precalculated ones"
# Reference values for G were taken from the spreadsheet Berman_Gibbs_Free_Energies.xlsx
#   (http://www.dewcommunity.org/uploads/4/1/7/6/41765907/sunday_afternoon_sessions__1_.zip accessed on 2017-10-03)
T <- c(100, 100, 1000, 1000)
P <- c(5000, 50000, 5000, 50000)

# anadalusite: an uncomplicated mineral (no transitions)
And_G <- c(-579368, -524987, -632421, -576834)
And <- subcrt("andalusite", T = T, P = P)$out[[1]]
expect_true(maxdiff(And$G, And_G) < 7.5, info = info)

# quartz: a mineral with polymorphic transitions
aQz_G <- c(-202800, -179757, -223864, -200109)
aQz <- subcrt("quartz", T = T, P = P)$out[[1]]
expect_true(maxdiff(aQz$G[-2], aQz_G[-2]) < 1.2, info = info)
# here, the high-P, low-T point suffers
expect_true(maxdiff(aQz$G[2], aQz_G[2]) < 1250, info = info)

# K-feldspar: this one has disordering effects
Kfs_G <- c(-888115, -776324, -988950, -874777)
Kfs <- subcrt("K-feldspar", T = T, P = P)$out[[1]]
expect_true(maxdiff(Kfs$G[1:2], Kfs_G[1:2]) < 5, info = info)
# we are less consistent with the reference values at high T
expect_true(maxdiff(Kfs$G[3:4], Kfs_G[3:4]) < 350, info = info)

info <- "Nonexistent or incomplete user data file is handled properly"
thermo("opt$Berman" = "XxXxXx.csv")
expect_error(berman("calcite"), "the file named in thermo\\(\\)\\$opt\\$Berman \\(XxXxXx.csv\\) does not exist", info = info)
thermo("opt$Berman" = system.file("extdata/Berman/testing/BA96_berman.csv", package = "CHNOSZ"))
expect_error(berman("xxx"), "Data for xxx not available. Please add it to", info = info)
thermo("opt$Berman" = NA)
expect_error(berman("xxx"), "Data for xxx not available. Please add it to your_data_file.csv", info = info)

info <- "NA values of P are handled"
sresult <- suppressWarnings(subcrt("H2O", T = seq(0, 500, 100)))
T <- sresult$out$water$T
P <- sresult$out$water$P
# This stopped with a error prior to version 1.1.3-37
bresult <- berman("quartz", T = convert(T, "K"), P = P)
expect_equal(sum(is.na(bresult$G)), 2, info = info)
# This also now works (producing the same NA values)
#subcrt("quartz", T = seq(0, 500, 100))

"NAs don't creep into calculations below 298.15 K for minerals with disorder parameters"
# 20191116
expect_false(any(is.na(subcrt("K-feldspar", P = 1, T = seq(273.15, 303.15, 5), convert = FALSE)$out[[1]]$G)), info = info)

# Get parameters for all available minerals
dat <- berman()
mineral <- unique(dat$name)

info <- "Properties of all minerals are computed without errors"
# Running this without error means that:
# - formulas for the minerals are found in thermo()$OBIGT
# - warning is produced for flourtremolite (GfPrTr(calc) >= 1000 J/mol different from GfPrTr(table))
# - use units = "cal" for comparison with Helgeson minerals below
expect_warning(properties <- lapply(mineral, berman, check.G = TRUE, units = "cal"),
               "fluortremolite", info = info)
# Save the results so we can use them in the next tests
Berman <- do.call(rbind, properties)

# Find the mineral data using Helgeson formulation
icr <- suppressMessages(info(mineral, "cr"))
add.OBIGT("SUPCRT92")
# NOTE: with check.it = TRUE (the default), this calculates Cp from the tabulated Maier-Kelley parameters
#Helgeson <- suppressMessages(info(icr, check.it = FALSE))
Helgeson <- suppressMessages(info(icr))

# Get the minerals that are present in *both* Berman and Helgeson versions
# All of these except rutile (Robie et al., 1978) reference Helgeson et al., 1978
iboth <- Helgeson$ref1 %in% c("HDNB78", "RHF78.4")
mineral <- mineral[iboth]
Berman <- Berman[iboth, ]
Helgeson <- Helgeson[iboth, ]

# Now we can compare Berman and Helgeson G, H, S, Cp, V
# Minerals with missing properties are not matched here
# (i.e. fluortremolite: no G and H in Helgeson data)

info <- "Berman and Helgeson tabulated properties have large differences for few minerals"
# Which minerals differ in DGf by more than 4 kcal/mol?
idiffG <- which(abs(Berman$G - Helgeson$G) > 4000)
DGf.list <- c("paragonite", "anthophyllite", "antigorite", "Ca-Al-pyroxene", "lawsonite", "margarite", "merwinite", "fluorphlogopite")
expect_true(identical(sort(mineral[idiffG]), sort(DGf.list)), info = info)

# Which minerals differ in DHf by more than 4 kcal/mol?
idiffH <- which(abs(Berman$H - Helgeson$H) > 4000)
DHf.list <- c("paragonite", "anthophyllite", "antigorite", "Ca-Al-pyroxene", "lawsonite", "margarite", "merwinite", "fluorphlogopite")
expect_true(identical(sort(mineral[idiffH]), sort(DHf.list)), info = info)

# Which minerals differ in S by more than 4 cal/K/mol?
idiffS <- which(abs(Berman$S - Helgeson$S) > 4)
DS.list <- c("epidote", "annite", "fluortremolite", "andradite")
expect_true(identical(sort(mineral[idiffS]), sort(DS.list)), info = info)

# Which minerals differ in Cp by more than 4 cal/K/mol?
idiffCp <- which(abs(Berman$Cp - Helgeson$Cp) > 4)
DCp.list <- c("glaucophane", "antigorite", "cristobalite,beta", "K-feldspar", "fluortremolite")
expect_true(identical(sort(mineral[idiffCp]), sort(DCp.list)), info = info)

# Which minerals differ in V by more than 1 cm^3/mol?
idiffV <- which(abs(Berman$V - Helgeson$V) > 1)
DV.list <- c("glaucophane", "anthophyllite", "antigorite", "chrysotile", "merwinite", "grunerite")
expect_true(identical(sort(mineral[idiffV]), sort(DV.list)), info = info)
