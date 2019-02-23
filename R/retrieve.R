# CHNOSZ/retrieve.R
# retrieve species with given elements
# 20190214 initial version
# 20190224 use ... for multiple arguments (define a chemical system)

retrieve <- function(..., include.electron = FALSE, include.groups = FALSE, state = NULL) {
  ## stoichiometric matrix
  # what are the formulas of species in the current database?
  formula <- thermo()$obigt$formula
  # get a previously calculated stoichiometric matrix, if it matches the current database
  stoich <- thermo()$stoich
  if(!is.null(stoich)) {
    # if it doesn't match the current database, don't use it
    if(!identical(rownames(stoich), formula)) stoich <- NULL
  }
  if(is.null(stoich)) {
    # Create the stoichiometric matrix for the current database
    # and suppress warning messages about missing elements
    message("retrieve: creating stoichiometric matrix")
    # NOTE: row names are the formulas, so we can detect if the database changes
    stoich <- suppressWarnings(i2A(formula))
    # store the stoichiometric matrix for later calculations
    thermo("stoich" = stoich)
  }

  ## species identification
  args <- list(...)
  ispecies <- numeric()
  for(elements in args) {
    not.present <- ! elements %in% colnames(stoich)
    if(any(not.present)) {
      if(sum(not.present)==1) stop(elements[not.present], " is not an element that is present in any species")
      else stop(paste(elements[not.present], collapse=", "), " are not elements that are present in any species")
    }
    # identify the species that have the elements
    has.elements <- rowSums(stoich[, elements, drop = FALSE] != 0) == length(elements)
    # which species are these (i.e. the species index)
    ispecies <- c(ispecies, which(has.elements))
    ispecies <- ispecies[!duplicated(ispecies)]
  }
  # for a chemical system, defined by multiple arguments, the species can not contain any _other_ elements
  if(length(args) > 1) {
    syselements <- unlist(args)
    isyselements <- colnames(thermo()$stoich) %in% syselements
    notsysstoich <- thermo()$stoich[, !isyselements]
    iother <- rowSums(notsysstoich[ispecies, ] != 0) > 0
    ispecies <- ispecies[!iother]
    # include the species for "Z" (charge)
    if(!include.electron) {
      ielectron <- names(ispecies) == "(Z-1)"
      ispecies <- ispecies[!ielectron]
    }
  }
  # exclude groups and filter states
  if(!include.groups) {
    igroup <- grepl("^\\[.*\\]$", thermo()$obigt$name[ispecies])
    ispecies <- ispecies[!igroup]
  }
  if(!is.null(state)) {
    istate <- thermo()$obigt$state[ispecies] %in% state
    ispecies <- ispecies[istate]
  }
  ispecies
}

