# CHNOSZ/duplex.R
# Combine diagrams for two metals
# 20200713 first version jmd

# Function to make a new "affinity" object from two diagrams;
# uses *secondary* balancing coefficients to combine the diagrams
duplex <- function(d1, d2, balance = NULL) {
  # Check that the basis species are the same
  if(!identical(d1$basis, d2$basis)) stop("basis species are not identical")
  # Check that the variables and their values are the same
  if(!identical(d1$vars, d2$vars)) stop("variable names are not identical")
  if(!identical(d1$vals, d2$vals)) stop("variable values are not identical")
  # Check that T and P are the same
  if(!identical(d1$T, d2$T)) stop("temperatures are not identical")
  if(!identical(d1$P, d2$P)) stop("pressures are not identical")
  # Check that we have plotvals and predominant (from diagram())
  if(is.null(d1$plotvals) | is.null(d1$predominant)) stop("d1 is missing 'plotvals' or 'predominant' components (not made by diagram()?)")
  if(is.null(d2$plotvals) | is.null(d2$predominant)) stop("d2 is missing 'plotvals' or 'predominant' components (not made by diagram()?)")

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

