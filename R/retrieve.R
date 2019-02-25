# CHNOSZ/retrieve.R
# retrieve species with given elements
# 20190214 initial version
# 20190224 use ... for multiple arguments (define a chemical system)

retrieve <- function(..., state = NULL, add.charge = TRUE, hide.groups = TRUE, req1 = FALSE) {
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
  # automatically add charge to a system 20190225
  if(add.charge & length(args) > 1) {
    if(!"Z" %in% unlist(args)) args <- c(args, "Z")
  }
  # for a numeric first argument, limit the result to only those species 20190225
  for(elements in args) {
    if(identical(elements, "all")) {
      ispecies <- 1:nrow(thermo()$obigt)
      names(ispecies) <- thermo()$obigt$formula
    } else {
      not.present <- ! elements %in% colnames(stoich)
      if(any(not.present)) {
        if(sum(not.present)==1) stop('"', elements[not.present], '" is not an element that is present in any species')
        else stop('"', paste(elements[not.present], collapse='", "'), '" are not elements that are present in any species')
      }
      # identify the species that have the elements
      has.elements <- rowSums(stoich[, elements, drop = FALSE] != 0) == length(elements)
      # which species are these (i.e. the species index)
      # for req1, remember the species containing the first element 20190225
      if(length(ispecies)==0) ispecies1 <- which(has.elements)
      ispecies <- c(ispecies, which(has.elements))
      ispecies <- ispecies[!duplicated(ispecies)]
    }
  }
  # for a chemical system, defined by multiple arguments, the species can not contain any _other_ elements
  if(length(args) > 1) {
    syselements <- unlist(args)
    isyselements <- colnames(thermo()$stoich) %in% syselements
    notsysstoich <- thermo()$stoich[, !isyselements]
    iother <- rowSums(notsysstoich[ispecies, ] != 0) > 0
    ispecies <- ispecies[!iother]
  }
  # keep only species that contain the first element
  if(req1) {
    ispecies <- intersect(ispecies1, ispecies)
    names(ispecies) <- thermo()$obigt$name[ispecies]
  }
  # exclude groups
  if(hide.groups) {
    igroup <- grepl("^\\[.*\\]$", thermo()$obigt$name[ispecies])
    ispecies <- ispecies[!igroup]
  }
  #if(hide.electron) {
  #  ielectron <- names(ispecies) == "(Z-1)"
  #  ispecies <- ispecies[!ielectron]
  #}
  #if(hide.proton) {
  #  iproton <- names(ispecies) == "H+"
  #  ispecies <- ispecies[!iproton]
  #}
  # filter states
  if(!is.null(state)) {
    istate <- thermo()$obigt$state[ispecies] %in% state
    ispecies <- ispecies[istate]
  }
  # for names, use e- instead of (Z-1)
  names(ispecies)[names(ispecies)=="(Z-1)"] <- "e-"
  ispecies
}
