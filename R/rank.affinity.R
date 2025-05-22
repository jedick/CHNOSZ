# Calculate average rank of affinities for species in different groups
# 20220416 jmd first version
# 20250522 rescale ranks

rank.affinity <- function(aout, groups, rescale = TRUE, percent = FALSE) {

  # Put the affinities into matrix form
  amat <- sapply(aout$values, as.numeric)
  # Calculate ranks
  # https://stackoverflow.com/questions/1412775/pmax-parallel-maximum-equivalent-for-rank-in-r
  arank <- apply(amat, 1, rank)

  # Count total number of species in all groups
  groups_vector <- unlist(groups)
  if(is.integer(groups_vector)) ntot <- length(groups_vector)
  if(is.numeric(groups_vector)) ntot <- length(groups_vector)
  if(is.logical(groups_vector)) ntot <- sum(groups_vector)

  # Get the bounds of average ranks for a group with one species
  min1 <- 1
  max1 <- ntot

  # Get the average rank for species in each group
  grank <- sapply(groups, function(group) {

    # Get number of species in this group
    if(inherits(group, "logical")) n <- sum(group)
    if(inherits(group, "integer")) n <- length(group)
    # Also handle indices classed as numeric 20250522
    if(inherits(group, "numeric")) n <- length(group)
    # Sum the ranks and divide by number of species
    rank_avg <- colSums(arank[group, , drop = FALSE]) / n

    if(rescale) {
      # Rescale ranks 20250522
      # Get the bounds of average ranks for a group with n species
      # Minimum is the average of 1..n for n species
      min <- sum(1:n) / n
      # The margin is the difference between the minimum and 1
      margin <- min - min1
      # Lower and upper bounds are symmetric, so we subtract the margin from total number of species to get the max
      max <- ntot - margin
      ## Factor to rescale average ranks from an n-species group to a 1-species group
      #scaling_factor <- (max1 - min1) / (max - min)
      ## To center the range, we have to subtract the margin on both sides
      #(rank_avg - 2 * margin) * scaling_factor

      # Build a linear model mapping from x (bounds of group with n species) to y (bounds of group with 1 species)
      x <- c(min, max)
      y <- c(min1, max1)
      rescale_lm <- lm(y ~ x)
      # Rescale average ranks with the linear model
      rank_avg <- predict(rescale_lm, data.frame(x = rank_avg))

    } else {
      rank_avg
    }

  })

  # Calculate average rank percentage 20240106
  if(percent) grank <- grank / rowSums(grank) * 100

  # Restore dims
  dims <- dim(aout$values[[1]])
  if(getRversion() < "4.1.0") {
    # Using 'simplify = FALSE' in R < 4.1.0 caused error: 3 arguments passed to 'dim<-' which requires 2
    glist <- lapply(lapply(apply(grank, 2, list), "[[", 1), "dim<-", dims)
  } else {
    # apply() got 'simplify' argument in R 4.1.0 20230313
    glist <- apply(grank, 2, "dim<-", dims, simplify = FALSE)
  }
  aout$values <- glist

  # Rename species to group names (for use by diagram())
  aout$species <- aout$species[1:length(groups), ]
  aout$species$name <- names(groups)
  # Label the object with our function name
  aout$fun <- "rank.affinity"
  aout

}
