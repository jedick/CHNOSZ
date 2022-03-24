# Load default settings for CHNOSZ
reset()

info <- "minimal usage of mod.OBIGT() creates usable data entries"
# We need at least a name and some property
expect_error(mod.OBIGT("test"), "species name and a property", info = info)
# A valid formula is needed
expect_warning(expect_error(mod.OBIGT("test", date = as.character(Sys.Date())), "is not a simple chemical formula", info = info),
             "please supply a valid chemical formula", info = info)
# The default state is aq
expect_message(itest <- mod.OBIGT("test", formula = "Z0", date = as.character(Sys.Date())), "added test\\(aq\\)", info = info)
# We should get NA values of G for a species with NA properties 
expect_true(all(is.na(subcrt(itest)$out[[1]]$G)), info = info)
# A single value of G comes through to subcrt
mod.OBIGT("test", G = 100)
expect_equal(subcrt("test", T = 25, P = 1)$out[[1]]$G, 100, info = info)
# Values for Cp and c1 integrate to the same values of G
G.Cp <- subcrt(mod.OBIGT(list(name = "test", S = 0, Cp = 100)))$out[[1]]$G
G.c1 <- subcrt(mod.OBIGT(list(name = "test", S = 0, c1 = 100)))$out[[1]]$G
expect_equal(G.Cp, G.c1, info = info)

info <- "mod.OBIGT() works with numeric argument (species index)"
# Test added 20200707
i1 <- info("CH4")
# Previously this failed with an error
i2 <- mod.OBIGT(i1, name = "methane")
expect_identical(i1, i2, info = info)
expect_identical(info(i2)$name, "methane", info = info)

# Test added 20220325
info <- "logK_to_OBIGT() errors for < 3 temperature values"
logB <- c(6.24, 6.02, 5.84, 5.97, 6.52)
species <- c("Co+2", "HS-", "CoHS+")
coeff <- c(-1, -1, 1)
T <- c(120, 150, 200, 250, 300)
P <- "Psat"
expect_error(logB_to_OBIGT(logB[1:2], species, coeff, T[1:2], P))
