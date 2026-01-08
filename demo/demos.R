# CHNOSZ/demo/demos.R
# Run all demos in the package

library(CHNOSZ)

demos <- c("references", "dehydration", "affinity", "NaCl", "density", 
  "ORP", "ionize", "buffer", "protbuff", "glycinate",
  "mosaic", "copper", "arsenic", "solubility", "gold", "contour", "sphalerite", "minsol",
  "Shh", "saturation", "adenine", "DEW", "lambda", "potassium", "TCA", "aluminum",
  "AD", "comproportionation", "Pourbaix", "E_coli", "yttrium", "rank.affinity", "uranyl",
  "sum_S", "MgATP", "rubisco_Zc", "phosphorylate")

for(i in 1:length(demos)) {
  # A message so the user knows where we are
  message("------------")
  message(paste("demos: running '", demos[i], "'", sep = ""))
  demo(demos[i], package = "CHNOSZ", character.only = TRUE, echo = FALSE, ask = FALSE)
}
