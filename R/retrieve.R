# CHNOSZ/retrieve.R
# retrieve species with given elements
# 20190214 initial version
# 20190224 define a chemical system using multiple arguments [defunct; use list() for chemical system]
# 20190304 update the stoichiometric matrix instead of doing a full recalculation when the database changes
# 20190305 use c() for combination of elements, list() for chemical system,
#          and add 'ligands' argument to retrieve element-bearing species

retrieve <- function(elements = NULL, ligands = NULL, state = NULL, T = NULL, P = "Psat", add.charge = TRUE, hide.groups = TRUE) {
  ## empty argument handling
  if(is.null(elements)) return(integer())

  ## stoichiometric matrix
  # what are the formulas of species in the current database?
  formula <- thermo()$OBIGT$formula
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

  ## generate error for missing element(s)
  allelements <- c(unlist(elements), unlist(ligands))
  not.present <- ! allelements %in% c(colnames(stoich), "all")
  if(any(not.present)) {
    if(sum(not.present)==1) stop('"', allelements[not.present], '" is not an element that is present in any species in the database')
    else stop('"', paste(allelements[not.present], collapse='", "'), '" are not elements that are present in any species in the database')
  }

  ## handle 'ligands' argument
  if(!is.null(ligands)) {

    # if 'ligands' is cr, liq, gas, or aq, use that as the state
    if(any(ligands %in% c("cr", "liq", "gas", "aq")) & length(ligands)==1) {
      state <- ligands
      ispecies <- retrieve(elements, add.charge = add.charge)
    } else {
      # include the element in the system defined by the ligands list
      ligands <- c(elements, as.list(ligands))
      # call retrieve() for each argument and take the intersection
      r1 <- retrieve(elements, add.charge = add.charge)
      r2 <- retrieve(ligands, add.charge = add.charge)
      ispecies <- intersect(r1, r2)
    }

  } else {

    ## species identification
    ispecies <- list()
    # automatically add charge to a system 20190225
    if(add.charge & is.list(elements) & !"Z" %in% elements) elements <- c(elements, "Z")
    # proceed element-by-element
    for(i in seq_along(elements)) {
      element <- unlist(elements[i])
      if(identical(element, "all")) {
        ispecies[[i]] <- 1:nrow(thermo()$OBIGT)
      } else {
        # identify the species that have the element
        has.element <- rowSums(stoich[, element, drop = FALSE] != 0) == 1
        ispecies[[i]] <- which(has.element)
      }
    }

    # now we have an list of ispecies (one vector for each element)
    # what we do next depends on whether the argument is a list() or c()
    if(is.list(elements)) {
      # for a chemical system, all species are included that do not contain any other elements
      ispecies <- unique(unlist(ispecies))
      ielements <- colnames(thermo()$stoich) %in% elements
      if(any(!ielements)) {
        otherstoich <- thermo()$stoich[, !ielements]
        iother <- rowSums(otherstoich[ispecies, ] != 0) > 0
        ispecies <- ispecies[!iother]
      }
    } else {
      # get species that have all the elements; the species must be present in each vector
      # Reduce() hint from https://stackoverflow.com/questions/27520310/union-of-intersecting-vectors-in-a-list-in-r
      ispecies <- Reduce(intersect, ispecies)
    }

  }

  # exclude groups
  if(hide.groups) {
    igroup <- grepl("^\\[.*\\]$", thermo()$OBIGT$name[ispecies])
    ispecies <- ispecies[!igroup]
  }
  # filter on state
  if(!is.null(state)) {
    istate <- thermo()$OBIGT$state[ispecies] %in% state
    ispecies <- ispecies[istate]
  }

  # require non-NA Delta G0 at specific temperature 20200825
  if(!is.null(T)) {
    G <- sapply(suppressMessages(subcrt(ispecies, T = T, P = P))$out, "[[", "G")
    ispecies <- ispecies[!is.na(G)]
  }

  # assign names; use e- instead of (Z-1)
  names(ispecies) <- thermo()$OBIGT$formula[ispecies]
  names(ispecies)[names(ispecies)=="(Z-1)"] <- "e-"
  # if there's nothing, don't give it a name
  if(length(ispecies)==0) ispecies <- integer()
  ispecies
}
