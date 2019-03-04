# CHNOSZ/retrieve.R
# retrieve species with given elements
# 20190214 initial version
# 20190224 use ... for multiple arguments (define a chemical system)
# 20190304 update the stoichiometric matrix instead of doing a full recalculation when the database changes

retrieve <- function(..., state = NULL, add.charge = TRUE, hide.groups = TRUE, req1 = FALSE) {
  ## stoichiometric matrix
  # what are the formulas of species in the current database?
  formula <- thermo()$obigt$formula
  # get a previously calculated stoichiometric matrix
  stoich <- thermo()$stoich
  # if it doesn't match the current database, update it
  if(!identical(rownames(stoich), formula)) {
    message("retrieve: updating stoichiometric matrix")
    # first put rows in the right order for the current database
    istoich <- match(formula, rownames(stoich))
    stoich <- stoich[na.omit(istoich), ]
    # deal with any missing formulas
    if(any(is.na(istoich))) {
      # convert the stoichiometric matrix to data frame with a 'formula' column
      oldstoich <- cbind(formula = rownames(stoich), as.data.frame(stoich))
      # get the stoichiometry of the additional formulas
      # and suppress warning messages about missing elements
      addstoich <- suppressWarnings(i2A(formula[is.na(istoich)]))
      addstoich <- cbind(formula = rownames(addstoich), as.data.frame(addstoich))
      # merge the old and added stoichiometries, and assign 0 for missing elements in either one
      newstoich <- merge(oldstoich, addstoich, all = TRUE)
      newstoich[is.na(newstoich)] <- 0
      # convert the data frame to matrix with rownames corresponding to formulas
      stoich <- as.matrix(newstoich[, 2:ncol(newstoich)])
      rownames(stoich) <- newstoich$formula
      # put rows in the right order for the current database
      istoich <- match(formula, rownames(stoich))
      stoich <- stoich[istoich, ]
    }
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
