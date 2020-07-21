# CHNOSZ/flatten.R
# Combine diagrams for two metals
# 20200713 first version jmd

# Function to combine two diagrams (simple overlay, no interaction) 20200717
# -- makes new "species" from all combinations of those in d1 and d2
flatten <- function(d1, d2) {
  check_d1_d2(d1, d2)

  # Index all combinations of species in d1 and d2
  i1 <- 1:nrow(d1$species)
  i2 <- 1:nrow(d2$species)
  combs <- expand.grid(i1, i2)

  # Get species rows for each combination
  s1 <- d1$species[combs[, 1], ]
  s2 <- d2$species[combs[, 2], ]
  # Make a new species data frame
  nbasis <- nrow(d1$basis)
  species <- s1[, 1:nbasis] + s2[, 1:nbasis]
  ispecies <- paste(s1$ispecies, s2$ispecies, sep = ",")
  logact <- paste(s1$logact, s2$logact, sep = ",")
  state <- paste(s1$state, s2$state, sep = ",")
  # Use names from diagram()
  if(is.expression(d1$names) & is.expression(d2$names)) {
    name <- lapply(1:nrow(combs), function(i) bquote(.(d1$names[[combs[i, 1]]])+.(d2$names[[combs[i, 2]]])))
    name <- unlist(lapply(name, deparse, width.cutoff = 500, control = NULL))
  } else if(is.expression(d1$names)) {
    name <- lapply(1:nrow(combs), function(i) bquote(.(d1$names[[combs[i, 1]]])+.(d2$names[combs[i, 2]])))
    name <- unlist(lapply(name, deparse, width.cutoff = 500, control = NULL))
  } else if(is.expression(d2$names)) {
    name <- lapply(1:nrow(combs), function(i) bquote(.(d1$names[combs[i, 1]])+.(d2$names[[combs[i, 2]]])))
    name <- unlist(lapply(name, deparse, width.cutoff = 500, control = NULL))
  } else name <- paste(d1$names[combs[, 1]], d2$names[combs[, 2]], sep="+")
  if(length(name) != nrow(combs)) stop("deparse()-ing expressions gives unequal length; try diagram(., format.names = FALSE)")
  species <- cbind(species, ispecies, logact, state, name)

  # Get affinities for each combination
  v1 <- d1$values[combs[, 1]]
  v2 <- d2$values[combs[, 2]]
  values <- Map("+", v1, v2)
  # Assign -Inf affinity where a species isn't predominant
  for(i in seq_along(values)) {
    i1 <- combs[i, 1]
    i2 <- combs[i, 2]
    ip1 <- d1$predominant == i1
    ip2 <- d2$predominant == i2
    ip12 <- ip1 & ip2
    values[[i]][!ip12] <- -Inf
  }

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
duplex <- function(d1, d2, balance = NULL) {
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

### unexported function ###

# Check that d1 and d2 can be combined
# Extracted from duplex() 20200717
check_d1_d2 <- function(d1, d2) {
  # Check that the basis species are the same
  if(!identical(d1$basis, d2$basis)) stop("basis species in objects 'd1' and 'd2' are not identical")
  # Check that the variables and their values are the same
  if(!identical(d1$vars, d2$vars)) stop("variable names in objects 'd1' and 'd2' are not identical")
  if(!identical(d1$vals, d2$vals)) stop("variable values in objects 'd1' and 'd2' are not identical")
  # Check that T and P are the same
  if(!identical(d1$T, d2$T)) stop("temperatures in objects 'd1' and 'd2' are not identical")
  if(!identical(d1$P, d2$P)) stop("pressures in objects 'd1' and 'd2' are not identical")
  # Check that we have plotvals and predominant (from diagram())
  if(is.null(d1$plotvals) | is.null(d1$predominant)) stop("object 'd1' is missing 'plotvals' or 'predominant' components (not made by diagram()?)")
  if(is.null(d2$plotvals) | is.null(d2$predominant)) stop("object 'd2' is missing 'plotvals' or 'predominant' components (not made by diagram()?)")
}
