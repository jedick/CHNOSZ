# CHNOSZ/demo/examples.R
# Run all examples in the package

library(CHNOSZ)

.ptime <- proc.time()
topics <- c("thermo", "examples",
  "util.array", "util.data", "util.expression", "util.legend", "util.plot",
  "util.formula", "util.misc", "util.seq", "util.units",
  "util.water", "taxonomy", "info", "retrieve", "add.OBIGT", "protein.info",
  "water", "IAPWS95", "subcrt", "Berman",
  "makeup", "basis", "swap.basis", "species", "affinity",
  "solubility", "equilibrate", 
  "diagram", "mosaic", "mix",
  "buffer", "nonideal", "NaCl",
  "add.protein", "ionize.aa", "EOSregress", "rank.affinity",
  "DEW", "logK.to.OBIGT", "stack_mosaic"
)
for(i in 1:length(topics)) example(topics[i], package = "CHNOSZ", character.only = TRUE, ask = FALSE)
cat("Time elapsed: ", proc.time() - .ptime, "\n")
