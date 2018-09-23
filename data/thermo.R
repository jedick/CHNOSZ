# CHNOSZ/data/thermo.R
# create or restore the 'thermo' data object

# create the CHNOSZ environment if it does not exist
if(!"CHNOSZ" %in% search()) {
  attach(NULL, name="CHNOSZ")
  message("data(thermo): attached environment \"CHNOSZ\"")
}

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
  # create thermo list
  thermo <- list(
    # as.is: keep character values as character and not factor
    opt = as.list(read.csv("opt.csv", as.is=TRUE)),
    element = read.csv("element.csv", as.is=1:3),
    #obigt = read.csv("OBIGT.csv", as.is=1:7),
    obigt = obigt,
    refs = read.csv("refs.csv", as.is=TRUE),
    buffers = read.csv("buffer.csv", as.is=1:3),
    protein = read.csv("protein.csv", as.is=1:4),
    groups = read.csv("groups.csv", row.names=1, check.names=FALSE),
    basis = NULL,
    species = NULL,
    opar = NULL
  )
  # place it in CHNOSZ environment
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
