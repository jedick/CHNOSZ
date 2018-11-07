# CHNOSZ/examples.R
# run examples from the help files, 
# and a function containing extra examples

examples <- function(save.png=FALSE) {
  # run all the examples in CHNOSZ documentation
  .ptime <- proc.time()
  topics <- c("thermo", "examples",
    "util.array", "util.blast", "util.data", "util.expression",
    "util.fasta", "util.formula", "util.matrix", "util.misc", "util.seq", "util.units",
    "util.water", "taxonomy", "info", "protein.info", "hkf", "water", "IAPWS95", "subcrt", "berman",
    "makeup", "basis", "swap.basis", "species", "affinity", "solubility", "equilibrate", 
    "diagram", "buffer", "nonideal", "NaCl", "add.protein", "protein", "ionize.aa", "yeast.aa",
    "objective", "revisit", "EOSregress", "wjd")
  plot.it <- FALSE
  if(is.character(save.png))
    png(paste(save.png,"%d.png",sep=""),width=500,height=500,pointsize=12)
  else if(save.png) plot.it <- TRUE
  for(i in 1:length(topics)) {
    if(plot.it) png(paste(topics[i],"%d.png",sep=""),width=500,height=500,pointsize=12)
    myargs <- list(topic=topics[i],ask=FALSE)
    do.call(example,myargs)
    if(plot.it) dev.off()
  }
  if(is.character(save.png)) dev.off()
  cat("Time elapsed: ", proc.time() - .ptime, "\n")
}

demos <- function(which=c("sources", "protein.equil", "affinity", "NaCl", "density", 
  "ORP", "revisit", "findit", "ionize", "buffer", "protbuff", "yeastgfp", "mosaic",
  "copper", "solubility", "gold", "wjd", "bugstab", "Shh", "activity_ratios",
  "adenine", "DEW", "lambda", "TCA", "go-IU", "bison"), save.png=FALSE) {
  # run one or more demos from CHNOSZ with ask=FALSE, and return the value of the last one
  for(i in 1:length(which)) {
    # say something so the user sees where we are
    message("------------")
    if(which[i]=="dehydration" & !save.png) {
      message("demos: skipping dehydration demo as save.png is FALSE")
      next 
    } else message(paste("demos: running '", which[i], "'", sep=""))
    if(save.png & !which[i]=="dehydration") {
      if(which[i]=="bugstab") png(paste(which[i], "%d.png", sep=""), width=700, height=500, pointsize=12)
      else png(paste(which[i], "%d.png", sep=""), width=500, height=500, pointsize=12)
    }
    out <- demo(which[i], package="CHNOSZ", character.only=TRUE, echo=FALSE, ask=FALSE)
    if(save.png & !which[i]=="dehydration") dev.off()
  }
  return(invisible(out))
}

