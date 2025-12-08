# CHNOSZ/examples.R
# Functions to run all examples and demos in the package

examples <- function(save.png = FALSE) {
  # Run examples in each of the CHNOSZ help pages
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
  plot.it <- FALSE
  if(is.character(save.png))
    png(paste(save.png, "%d.png", sep = ""), width = 500, height = 500, pointsize = 12)
  else if(save.png) plot.it <- TRUE
  for(i in 1:length(topics)) {
    if(plot.it) png(paste(topics[i], "%d.png", sep = ""), width = 500, height = 500, pointsize = 12)
    myargs <- list(topic = topics[i], ask = FALSE)
    do.call(example, myargs)
    if(plot.it) dev.off()
  }
  if(is.character(save.png)) dev.off()
  cat("Time elapsed: ", proc.time() - .ptime, "\n")
}

demos <- function(which = c("references", "dehydration", "affinity", "NaCl", "density", 
  "ORP", "ionize", "buffer", "protbuff", "glycinate",
  "mosaic", "copper", "arsenic", "solubility", "gold", "contour", "sphalerite", "minsol",
  "Shh", "saturation", "adenine", "DEW", "lambda", "potassium", "TCA", "aluminum",
  "AD", "comproportionation", "Pourbaix", "E_coli", "yttrium", "rank.affinity", "uranyl",
  "sum_S", "MgATP", "rubisco_Zc", "phosphorylate"),
  save.png = FALSE) {
  # Run one or more demos from CHNOSZ with ask = FALSE, and return the value of the last one
  out <- NULL
  for(i in 1:length(which)) {
    # A message so the user knows where we are
    message("------------")
    message(paste("demos: running '", which[i], "'", sep = ""))
    width <- 500
    height <- 500
    pointsize <- 12
    if(which[i] == "comproportionation") width <- 600
    if(which[i] == "phosphorylate") {width <- 800; height <- 600; pointsize <- 16}
    if(save.png) png(paste(which[i], "%d.png", sep = ""), width = width, height = height, pointsize = pointsize)
    out <- demo(which[i], package = "CHNOSZ", character.only = TRUE, echo = FALSE, ask = FALSE)
    if(save.png) dev.off()
  }
  return(invisible(out))
}

