# CHNOSZ/mix.R
# Combine diagrams for two metals
# 20200713 first version jmd

## If this file is interactively sourced, the following are also needed to provide unexported functions:
#source("equilibrate.R")

# Function to combine two diagrams (simple overlay, no interaction) 20200717
# -- makes new "species" from all combinations of those in d1 and d2
mash <- function(d1, d2) {
  # It's just mixing d1 and d2 (two single-metal diagrams) without adding d3 (bimetal) 20200722
  mix(d1, d2)
}

# Mix the systems to include bimetal species 20200721
# 'parts' gives number of moles of each metal in 1 mole of the mixture
mix <- function(d1, d2, d3 = NULL, parts = c(1, 1), .balance = NULL) {
  if(is.null(.balance)) check_d1_d2(d1, d2)
  else check_d1_d3(d1, d3)

  # Handle mixing here
  if(!is.null(d3) & is.null(.balance)) {
    # Mix d1 (e.g. Fe only) and d2 (e.g. V only)
    m12 <- mix(d1, d2, parts = parts)
    # Mix d1 (Fe) and d3 (FeV bimetallic species)
    # - '.balance' is defined to add d3 to get required amount of V
    # - second 'd3' in argument list is used only for check_d1_d3()
    m13 <- mix(d1, d3, d3, parts = parts, .balance = d2$balance)
    # Mix d2 (V) and d3 (FeV)
    # - '.balance' is defined to add d3 to get required amount of Fe
    # - d2 (V) is first, so we need to reverse the 'parts' values
    m23 <- mix(d2, d3, d3, parts = rev(parts), .balance = d1$balance)
    # Merge all the species and affinity values
    species <- rbind(m12$species, m13$species, m23$species)
    values <- c(m12$values, m13$values, m23$values)
    if(nrow(d3$species) > 1) {
      # Mix d3 with itself (combinations of bimetallic species)
      m33 <- mix(d3, d3, d3, parts = parts, .balance = c(d1$balance, d2$balance))
      # Merge all the species and affinity values
      species <- rbind(species, m33$species)
      values <- c(values, m33$values)
    } 
    # Remove duplicates
    # (i.e. bimetallic species that exactly match the composition in 'parts'
    # and therefore appear in multiple combinations with
    # mono-metallic species that have zero mole fractions)
    idup <- duplicated(species$name)
    species <- species[!idup, ]
    values <- values[!idup]
    # Put together the output
    anew <- m12
    anew$species <- species
    anew$values <- values
    return(anew)
  }

  # Index all combinations of species in d1 and d2
  i1 <- 1:nrow(d1$species)
  i2 <- 1:nrow(d2$species)
  combs <- expand.grid(i1, i2)
  # For bimetallic species (m33)
  if(length(.balance)==2) {
    # Remove combinations of a species with itself
    combs <- combs[combs[, 1] != combs[, 2], ]
    # Remove duplicated combinations
    isdup <- duplicated(paste(apply(combs, 1, sort)[1, ], apply(combs, 1, sort)[2, ]))
    combs <- combs[!isdup, ]
  }
  # Get species rows for each combination
  s1 <- d1$species[combs[, 1], ]
  s2 <- d2$species[combs[, 2], ]
  # Get balancing coefficients for each combination
  b1 <- d1$n.balance[combs[, 1]]
  b2 <- d2$n.balance[combs[, 2]]
  # Locate the columns for the two metals in the formation reactions
  ibal <- match(c(d1$balance, d2$balance), colnames(s1))
  # With non-null '.balance', d2 is really d3 so we calculate the needed balancing coefficients
  if(!is.null(.balance)) {
    if(length(.balance)==2) {
      # For combinations of bimetallic species (m33)
      b1 <- suppressMessages(balance(d1, .balance[1])$n.balance[combs[, 1]])
      b2 <- suppressMessages(balance(d2, .balance[2])$n.balance[combs[, 2]])
      ibal <- match(.balance, colnames(s1))
    } else {
      b2 <- suppressMessages(balance(d2, .balance)$n.balance[combs[, 2]])
      ibal <- match(c(d1$balance, .balance), colnames(s1))
    }
  }
  # Solve for the mole fractions of each species that give the required mixture
  isingular <- frac1 <- frac2 <- numeric()
  for(i in 1:nrow(combs)) {
    # The stoichiometric matrix
    A <- matrix(as.numeric(c(s1[i, ibal], s2[i, ibal])), nrow = 2, byrow = TRUE)
    x <- tryCatch(solve(t(A), parts), error = function(e) e)
    # Check if we have a singular combination (e.g. FeV and FeVO4)
    if(inherits(x, "error")) {
      isingular <- c(isingular, i)
    } else {
      frac1 <- c(frac1, x[1])
      frac2 <- c(frac2, x[2])
    }
  }
  if(length(isingular) > 0) {
    if(length(isingular) > 1) stxt <- "s" else stxt <- ""
    message(paste0("mix: removing ", length(isingular), " combination", stxt, " with a singular stoichiometric matrix"))
    combs <- combs[-isingular, ]
    s1 <- s1[-isingular, ]
    s2 <- s2[-isingular, ]
  }
  # Note that some of frac1 or frac2 might be < 0 ... we remove these combinations below

  # Make a new species data frame starting with sum of scaled formation reactions
  nbasis <- nrow(d1$basis)
  species <- s1[, 1:nbasis] * frac1 + s2[, 1:nbasis] * frac2
  # Concatenate ispecies, logact, state
  ispecies <- paste(s1$ispecies, s2$ispecies, sep = ",")
  logact <- paste(s1$logact, s2$logact, sep = ",")
  state <- paste(s1$state, s2$state, sep = ",")
  # Omit species with zero mole fraction from concatenated values 20200722
  ispecies[frac1 == 0] <- s2$ispecies[frac1 == 0]
  logact[frac1 == 0] <- s2$logact[frac1 == 0]
  state[frac1 == 0] <- s2$state[frac1 == 0]
  ispecies[frac2 == 0] <- s1$ispecies[frac2 == 0]
  logact[frac2 == 0] <- s1$logact[frac2 == 0]
  state[frac2 == 0] <- s1$state[frac2 == 0]
  # Use names from diagram()
  if(is.expression(d1$names) | is.expression(d2$names)) {
    # Convert non-expressions to lists so we can use [[ indexing below
    if(!is.expression(d1$names)) d1$names <- as.list(d1$names)
    if(!is.expression(d2$names)) d2$names <- as.list(d2$names)
    name <- lapply(1:nrow(combs), function(i) {
      # Don't include names of species that are not present (zero mole fractions)
      if(frac1[i] == 0 & frac2[i] == 0) bquote("")
      else if(frac1[i] == 0) bquote(.(d2$names[[combs[i, 2]]]))
      else if(frac2[i] == 0) bquote(.(d1$names[[combs[i, 1]]]))
      else {
        # Put the names together, with the species from d2 first if it is has a higher mole fraction 20200722
        if(frac2[i] > frac1[i]) bquote(.(d2$names[[combs[i, 2]]])+.(d1$names[[combs[i, 1]]]))
        else bquote(.(d1$names[[combs[i, 1]]])+.(d2$names[[combs[i, 2]]]))
      }
    })
    name <- unlist(lapply(name, deparse, width.cutoff = 500, control = NULL))
    if(length(name) != nrow(combs)) stop("deparse()-ing expressions gives unequal length; try diagram(., format.names = FALSE)")
  } else {
    # Plain text names
    name <- sapply(1:nrow(combs), function(i) {
      if(frac1[i] <= 0) d2$names[combs[i, 2]]
      else if(frac2[i] <= 0) d1$names[combs[i, 1]]
      else paste(d1$names[combs[i, 1]], d2$names[combs[i, 2]], sep="+")
    })
  }
  species <- cbind(species, ispecies, logact, state, name, stringsAsFactors = FALSE)

  # Get affinities for each combination of species
  v1 <- d1$values[combs[, 1]]
  v2 <- d2$values[combs[, 2]]
  # Scale affinities by mole fractions computed for the mixture
  v1 <- Map("*", v1, as.list(frac1))
  v2 <- Map("*", v2, as.list(frac2))
  # Add together the scaled affinities
  values <- Map("+", v1, v2)
  ipredominant <- logical(length(values))
  if(length(.balance)==2) {
    # For combinations of bimetallic species (m33), don't do predominance masking
    # (there are no mono-metallic species to look at)
    ipredominant[] <- TRUE
  } else {
    # Loop over combinations to find predominant species in the single-metal diagram(s)
    for(i in seq_along(values)) {
      # Get predominant species in first diagram
      i1 <- combs[i, 1]
      ip1 <- d1$predominant == i1
      # If the mole fraction is zero, it is predominant by definition
      # (this allows a single bimetallic species (that is paired with
      #  the zero-mole-fraction mono-metallic species) to be formed)
      if(frac1[i] == 0) ip1[] <- TRUE
      if(is.null(.balance)) {
        # Get predominant species in second diagram
        i2 <- combs[i, 2]
        ip2 <- d2$predominant == i2
        # If the mole fraction is zero, it is predominant by definition
        # (this allows us to recover the first single-element diagram with frac = c(1, 0))
        if(frac2[i] == 0) ip2[] <- TRUE
        # Assign -Inf affinity where any species isn't predominant in the corresponding single-metal diagram
        ip12 <- ip1 & ip2
        values[[i]][!ip12] <- -Inf
        if(any(ip12)) ipredominant[i] <- TRUE
      } else {
        # For non-null '.balance', only the first diagram is for a single metal
        values[[i]][!ip1] <- -Inf
        if(any(ip1)) ipredominant[i] <- TRUE
      }
    }
  }
  # Remove combinations that:
  # 1) involve a mono-metallic species with no predominance field in the corresponding single-metal diagram or
  # 2) have a negative mole fraction of any species
  inotnegative <- frac1 >= 0 & frac2 >= 0
  values <- values[ipredominant & inotnegative]
  species <- species[ipredominant & inotnegative, ]

  # Use d1 as a template for the new affinity object
  anew <- d1[1:11]
  # Insert combined results
  anew$species <- species
  anew$values <- values
  # We don't have sout (results from subcrt()) for the combined "species"
  anew$sout <- NULL
  anew
}

# Function to make a new "affinity" object from two diagrams 20200713
# -- uses *secondary* balancing coefficients to combine the diagrams
rebalance <- function(d1, d2, balance = NULL) {
  check_d1_d2(d1, d2)

  # Combine the species data frames
  species <- rbind(d1$species, d2$species)
  # Combine the sout objects (results from subcrt())
  only2 <- !d2$sout$species$ispecies %in% d1$sout$species$ispecies
  sout <- d1$sout
  sout$species <- rbind(sout$species, d2$sout$species[only2, ])
  sout$out <- c(sout$out, d2$sout$out[only2])
  # Combine the affinity values divided by the *primary*
  # balancing coefficients ("plotvals" from diagram())
  values <- c(d1$plotvals, d2$plotvals)

  # Use d1 as a template for the new affinity object
  anew <- d1[1:11]
  # Insert combined results
  anew$species <- species
  anew$sout <- sout
  anew$values <- values

  # Figure out the *secondary* balancing coefficients
  n.balance <- balance(anew, balance = balance)$n.balance
  # In the Fe-Cu-S-O-H example all the coefficients on H+ are negative
  if(all(n.balance < 0)) n.balance <- -n.balance
  n1 <- nrow(d1$species)
  n.balance.1 <- n.balance[1:n1]
  n.balance.2 <- n.balance[(n1+1):length(n.balance)]

  # Make empty matrices to hold affinities and balancing coefficients
  a1 <- d1$values[[1]]
  a1[] <- NA
  b2 <- a2 <- b1 <- a1
  # Get the affinities (per mole of species, not divided by any balancing coefficients)
  # and the secondary balancing coefficients for the predominant species in each diagram
  p1 <- d1$predominant
  for(ip in unique(as.vector(p1))) {
    a1[p1 == ip] <- d1$values[[ip]][p1 == ip]
    b1[p1 == ip] <- n.balance.1[ip]
  }
  p2 <- d2$predominant
  for(ip in unique(as.vector(p2))) {
    a2[p2 == ip] <- d2$values[[ip]][p2 == ip]
    b2[p2 == ip] <- n.balance.2[ip]
  }
  # Divide the affinities by the secondary balancing coefficients
  ab1 <- a1 / b1
  ab2 <- a2 / b2
  # Identify the species with the highest affinity (predominant in the *secondary* reactions)
  i1 <- ab1 > ab2
  # Suppress non-predominant species at each grid point
  for(i in 1:n1) anew$values[[i]][!i1] <- -Inf
  for(i in (n1+1):length(n.balance)) anew$values[[i]][i1] <- -Inf

  anew

}


### Unexported functions ###

# Check that d1 and d2 can be combined
# Extracted from duplex() (now rebalance()) 20200717
check_d1_d2 <- function(d1, d2) {
  # Check that the basis species are the same
  if(!identical(d1$basis, d2$basis)) stop(paste0("basis species in objects 'd1' and 'd2' are not identical"))
  # Check that the variables and their values are the same
  if(!identical(d1$vars, d2$vars)) stop(paste0("plot variables in objects 'd1' and 'd2' are not identical"))
  if(!identical(d1$vals, d2$vals)) stop(paste0("plot ranges in objects 'd1' and 'd2' are not identical"))
  # Check that T and P are the same
  if(!identical(d1$T, d2$T)) stop(paste0("temperatures in objects 'd1' and 'd2' are not identical"))
  if(!identical(d1$P, d2$P)) stop(paste0("pressures in objects 'd1' and 'd2' are not identical"))
  # Check that we have plotvals and predominant (from diagram())
  if(is.null(d1$plotvals) | is.null(d1$predominant)) stop("object 'd1' is missing 'plotvals' or 'predominant' components (not made by diagram()?)")
  if(is.null(d2$plotvals) | is.null(d2$predominant)) stop(paste0("object 'd2' is missing 'plotvals' or 'predominant' components (not made by diagram()?)"))
}

# Check that d1 and d3 can be combined
check_d1_d3 <- function(d1, d3) {
  # Check that the basis species are the same
  if(!identical(d1$basis, d3$basis)) stop(paste0("basis species in objects 'd1' and 'd3' are not identical"))
  # Check that the variables and their values are the same
  if(!identical(d1$vars, d3$vars)) stop(paste0("plot variables in objects 'd1' and 'd3' are not identical"))
  if(!identical(d1$vals, d3$vals)) stop(paste0("plot ranges in objects 'd1' and 'd3' are not identical"))
  # Check that T and P are the same
  if(!identical(d1$T, d3$T)) stop(paste0("temperatures in objects 'd1' and 'd3' are not identical"))
  if(!identical(d1$P, d3$P)) stop(paste0("pressures in objects 'd1' and 'd3' are not identical"))
  # Check that we have plotvals and predominant (from diagram())
  if(is.null(d1$plotvals) | is.null(d1$predominant)) stop("object 'd1' is missing 'plotvals' or 'predominant' components (not made by diagram()?)")
  if(is.null(d3$plotvals) | is.null(d3$predominant)) stop(paste0("object 'd3' is missing 'plotvals' or 'predominant' components (not made by diagram()?)"))
}
