# CHNOSZ/stack_mosaic.R

# Function to create mosaic stack 20220617
# Adapted from vignettes/multi-metal.Rmd
# col: Colors for species1, species2, and species12
#   (default of NA for col[3] means to plot species12 boundaries with color for species2)
# ...: Arguments for mosaic() (including affinity() arguments)
stack_mosaic <- function(bases, species1, species2, species12, names = NULL, col = list(4, 3, 6), col.names = list(4, 3, 6),
  fill = NULL, dx = list(0, 0, 0), dy = list(0, 0, 0), srt = list(0, 0, 0), lwd = list(1, 1, 1), lty = list(1, 1, 1),
  loga_aq = NULL, plot.it = TRUE, ...) {

  # Default is to use semi-transparent fill for bimetallic species
  if(is.null(fill)) fill <- list(NA, NA, add.alpha(col.names[3], "50"))

  # Load species1 (first metal-bearing species)
  isp1 <- species(species1)
  # Set activities 20220722
  if(!is.null(loga_aq)) {
    iaq1 <- which(isp1$state == "aq")
    if(length(iaq1) > 0) species(iaq1, loga_aq)
    loga_aq_arg <- c(NA, loga_aq)
  } else loga_aq_arg <- NULL
  # Calculate affinity of species1 while speciating bases (e.g. aqueous S species)
  mosaic1 <- mosaic(bases, loga_aq = loga_aq_arg, ...)
  # Show predominance fields
  diagram1 <- diagram(mosaic1$A.species, names = names[[1]], col = col[[1]], col.names = col.names[[1]], fill = fill[[1]],
    dx = dx[[1]], dy = dy[[1]], srt = srt[[1]], lwd = lwd[[1]], lty = lty[[1]], plot.it = plot.it)

  # Load species12 (bimetallic species) and species2 (second metal-bearing species)
  isp2 <- species(c(species12, species2))
  # Set activities 20220722
  if(!is.null(loga_aq)) {
    iaq2 <- which(isp2$state == "aq")
    if(length(iaq2) > 0) species(iaq2, loga_aq)
  }
  # Speciate bases again (NULL)
  # Take the predominant members of species1 (diagram1$predominant)
  mosaic2 <- mosaic(list(bases, species1), stable = list(NULL, diagram1$predominant), loga_aq = loga_aq_arg, ...)

  # Set colors
  col <- c(rep_len(col[[3]], length(species12)), rep_len(col[[2]], length(species2)))
  col.names <- c(rep_len(col.names[[3]], length(species12)), rep_len(col.names[[2]], length(species2)))
  # For NULL names, use the species names
  if(is.null(names[[3]])) names[[3]] <- species12
  if(is.null(names[[2]])) names[[2]] <- species2
  names <- c(names[[3]], names[[2]])
  # Set other parameters
  dx <- c(rep_len(dx[[3]], length(species12)), rep_len(dx[[2]], length(species2)))
  dy <- c(rep_len(dy[[3]], length(species12)), rep_len(dy[[2]], length(species2)))
  srt <- c(rep_len(srt[[3]], length(species12)), rep_len(srt[[2]], length(species2)))
  lwd <- c(rep_len(lwd[[3]], length(species12)), rep_len(lwd[[2]], length(species2)))
  lty <- c(rep_len(lty[[3]], length(species12)), rep_len(lty[[2]], length(species2)))
  fill <- c(rep_len(fill[[3]], length(species12)), rep_len(fill[[2]], length(species2)))
  diagram2 <- diagram(mosaic2$A.species, add = TRUE, names = names, col = col, col.names = col.names, fill = fill,
    dx = dx, dy = dy, srt = srt, lwd = lwd, lty = lty, plot.it = plot.it)

  out <- list(diagram1, diagram2)
  invisible(out)

}
