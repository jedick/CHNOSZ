# CHNOSZ/demo/demos_png.R
# Run all demos in the package, saving output to PNG files

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
  width <- 500
  height <- 500
  pointsize <- 12
  if(demos[i] == "comproportionation") width <- 600
  if(demos[i] == "phosphorylate") {width <- 800; height <- 600; pointsize <- 16}
  png(paste(demos[i], "%d.png", sep = ""), width = width, height = height, pointsize = pointsize)
  demo(demos[i], package = "CHNOSZ", character.only = TRUE, echo = FALSE, ask = FALSE)
  dev.off()
}
