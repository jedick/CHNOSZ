# CHNOSZ/protein.info.R

# Calculate formulas and summarize properties of proteins
# pinfo: find rownumber in thermo()$protein
# protein.length: lengths of the indicated proteins
# protein.formula: chemical makeup of the indicated proteins
# protein.OBIGT: perform group additivity calculations
# protein.basis: coefficients of basis species in formation reactions of [ionized] proteins [residues]

pinfo <- function(protein, organism = NULL, residue = FALSE, regexp = FALSE) {
  # Return the `protein` (possibly per residue) for:
  #   dataframe `protein`
  # Return the rownumber(s) of thermo()$protein for:
  #   character `protein`, e.g. LYSC_CHICK
  #   character `protein` and `organism`, e.g. 'LYSC', 'CHICK'
  # Return the row(s) of thermo()$protein (possibly per residue) for:
  #   numeric `protein` (the rownumber itself)
  if(is.data.frame(protein)) out <- protein
  if(is.numeric(protein)) {
    t_p <- get("thermo", CHNOSZ)$protein
    # Drop NA matches to thermo()$protein
    iproteins <- 1:nrow(t_p)
    protein[!protein %in% iproteins] <- NA
    # Get amino acid counts
    out <- t_p[protein, ]
  }
  if(is.data.frame(protein) | is.numeric(protein)) {
    # Compute per-residue counts if requested
    if(residue) out[, 5:25] <- out[, 5:25]/rowSums(out[, 6:25])
  } else {
    t_p <- get("thermo", CHNOSZ)$protein
    # Search for protein by regular expression
    if(regexp) {
      iprotein <- grepl(protein, t_p$protein)
      iorganism <- iprotein
      if(!is.null(organism)) iorganism <- grepl(organism, t_p$organism)
      iprotein <- which(iprotein & iorganism)
    } else {
      # Search for protein or protein_organism in thermo()$protein
      t_p_names <- paste(t_p$protein, t_p$organism, sep = "_")
      if(is.null(organism)) my_names <- protein
      else my_names <- paste(protein, organism, sep = "_")
      iprotein <- match(my_names, t_p_names)
    }
    out <- iprotein
  }
  out
}

protein.formula <- function(protein, organism = NULL, residue = FALSE) {
  # Return a matrix with chemical formulas of proteins
  aa <- pinfo(pinfo(protein, organism))
  rf <- group.formulas()
  out <- as.matrix(aa[, 5:25]) %*% as.matrix(rf)
  if(residue) out <- out / rowSums(aa[, 6:25])
  row.names(out) <- make.unique(paste(aa$protein, aa$organism, sep = "_"))
  return(out)
}

protein.length <- function(protein, organism = NULL) {
  # Calculate the length(s) of proteins
  aa <- pinfo(pinfo(protein, organism))
  # Use rowSums on the columns containing amino acid counts
  pl <- as.numeric(rowSums(aa[, 6:25]))
  return(pl)
}

protein.OBIGT <- function(protein, organism = NULL, state = thermo()$opt$state) {
  # Display and return the properties of
  # proteins calculated from amino acid composition
  aa <- pinfo(pinfo(protein, organism))
  # The names of the protein backbone groups depend on the state
  # [UPBB] for aq or [PBB] for cr
  if(state == "aq") bbgroup <- "UPBB" else bbgroup <- "PBB"
  # Names of the AABB, sidechain and protein backbone groups
  groups <- c("AABB", colnames(aa)[6:25], bbgroup)
  # Put brackets around the group names
  groups <- paste("[", groups, "]", sep = "")
  # The rownumbers of the groups in thermo()$OBIGT
  groups_state <- paste(groups, state)
  OBIGT <- get("thermo", CHNOSZ)$OBIGT
  OBIGT_state <- paste(OBIGT$name, OBIGT$state)
  igroup <- match(groups_state, OBIGT_state)
  # The properties are in columns 10-22 of thermo()$OBIGT
  groupprops <- OBIGT[igroup, 10:22]
  # The elements in each of the groups
  groupelements <- i2A(igroup)
  # A function to work on a single row of aa
  eosfun <- function(aa) {
    # Numbers of groups: chains [=AABB], sidechains, protein backbone
    nchains <- as.numeric(aa[, 5])
    length <- sum(as.numeric(aa[, 6:25]))
    npbb <- length - nchains
    ngroups <- c(nchains, as.numeric(aa[, 6:25]), npbb)
    # The actual adding and multiplying of thermodynamic properties
    # Hmm. seems like we have to split up the multiplication/transposition
    # operations to get the result into multiple columns. 20071213
    eos <- t(data.frame(colSums(groupprops * ngroups)))
    # To get the formula, add up and round the group compositions 20090331
    f.in <- round(colSums(groupelements * ngroups), 3)
    # Take out any elements that don't appear (sometimes S)
    f.in <- f.in[f.in != 0]
    # Turn it into a formula
    f <- as.chemical.formula(f.in)
    # Now the species name
    name <- paste(aa$protein, aa$organism, sep = "_")
    # Tell the user about it
    message("protein.OBIGT: found ", appendLF = FALSE)
    message(name, " (", f, ", ", appendLF = FALSE)
    message(round(length, 3), " residues)")
    ref <- aa$ref
    # Include 'model' column 20220919
    model <- ifelse(state == "aq", "HKF", "CGL")
    header <- data.frame(name = name, abbrv = NA, formula = f, state = state, ref1 = ref, ref2 = NA,
      date = NA, model = model, E_units = "cal", stringsAsFactors = FALSE)
    eosout <- cbind(header, eos)
    return(eosout)
  }
  # Loop over each row of aa
  out <- lapply(1:nrow(aa), function(i) eosfun(aa[i, ]))
  out <- do.call(rbind, out)
  rownames(out) <- NULL
  return(out)
}

protein.basis <- function(protein, T = 25, normalize = FALSE) {
  # 20090902 Calculate the coefficients of basis species in reactions
  # to form proteins (possibly per normalized by length) listed in protein
  # 20120528 Renamed protein.basis from residue.info
  # What are the elemental compositions of the proteins
  aa <- pinfo(pinfo(protein))
  pf <- protein.formula(aa)
  # What are the coefficients of the basis species in the formation reactions
  sb <- species.basis(pf)
  # Calculate ionization states if H+ is a basis species
  thermo <- get("thermo", CHNOSZ)
  iHplus <- match("H+", rownames(thermo$basis))
  if(!is.na(iHplus)) {
    pH <- -thermo$basis$logact[iHplus]
    Z <- ionize.aa(aa, T = T, pH = pH)[1, ]
    sb[, iHplus] <- sb[, iHplus] + Z
  }
  # Compute per length-normalized coefficients if requested
  if(normalize) {
    # get lengths of proteins
    plen <- protein.length(aa)
    sb <- sb/plen
  }
  # Return the result
  return(sb)
}

