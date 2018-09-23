# CHNOSZ/data/OBIGT.R
# default thermodynamic data (OBIGT) in thermo

# we only work if the CHNOSZ environment exists
if(!"CHNOSZ" %in% search()) {
  message("data(OBIGT): please run data(thermo) first")
} else {

  local({
    # create obigt data frame
    sources_aq <- paste0(c("H2O", "inorganic", "organic", "biotic"), "_aq")
    sources_cr <- paste0(c("inorganic", "organic", "Berman"), "_cr")
    sources_liq <- paste0(c("organic"), "_liq")
    sources_gas <- paste0(c("inorganic", "organic"), "_gas")
    sources <- c(sources_aq, sources_cr, sources_liq, sources_gas)
    OBIGTdir <- system.file("extdata/OBIGT/", package="CHNOSZ")
    # need explicit "/" for Windows
    sourcefiles <- paste0(OBIGTdir, "/", c(sources_aq, sources_cr, sources_liq, sources_gas), ".csv")
    sourcefiles[!sources=="Berman_cr"] <- paste0(sourcefiles[!sources=="Berman_cr"], ".xz")
    datalist <- lapply(sourcefiles, read.csv, as.is=TRUE)
    obigt <- do.call(rbind, datalist)
    # get thermo from CHNOSZ environment
    thermo <- get("thermo", "CHNOSZ")
    # set obigt component
    thermo$obigt <- obigt
    # place thermo in CHNOSZ environment
    assign("thermo", thermo, "CHNOSZ")
  })

  # give a summary of some of the data
  message(paste("thermo$obigt:",
    nrow(thermo$obigt[thermo$obigt$state=="aq",]),
    "aqueous,", nrow(thermo$obigt), "total species"))

  # warn if there are duplicated species
  local({
    idup <- duplicated(paste(thermo$obigt$name, thermo$obigt$state))
    if(any(idup)) warning("thermo$obigt: duplicated species: ", 
      paste(thermo$obigt$name[idup], "(", thermo$obigt$state[idup], ")", sep="", collapse=" "))
  })

}
