# Calculate normalized sum of ranking of affinities for species in designated groups
# 20220416 jmd first version

rank.affinity <- function(aout, groups, percent = TRUE) {
  # Put the affinities into matrix form
  amat <- sapply(aout$values, as.numeric)
  # Calculate ranks
  # https://stackoverflow.com/questions/1412775/pmax-parallel-maximum-equivalent-for-rank-in-r
  arank <- apply(amat, 1, rank)
  # Get the normalized ranks for each group
  grank <- sapply(groups, function(group) {
    # Sum the ranks for this group and divide by number of species in the group
    if(inherits(group, "logical")) n <- sum(group)
    if(inherits(group, "integer")) n <- length(group)
    colSums(arank[group, ]) / n
  })
  # Calculate rank-sum percentage 20240106
  if(percent) grank <- grank / rowSums(grank) * 100
  # Restore dims
  dims <- dim(aout$values[[1]])
  # apply() got 'simplify' argument in R 4.1.0 20230313
  # Using 'simplify = FALSE' in R < 4.1.0 caused error: 3 arguments passed to 'dim<-' which requires 2
  if(getRversion() < "4.1.0") glist <- lapply(lapply(apply(grank, 2, list), "[[", 1), "dim<-", dims)
  else glist <- apply(grank, 2, "dim<-", dims, simplify = FALSE)
  aout$values <- glist
  # Rename species to group names (for use by diagram())
  aout$species <- aout$species[1:length(groups), ]
  aout$species$name <- names(groups)
  # "Sign" the object with our function name
  aout$fun <- "rank.affinity"
  aout
}
