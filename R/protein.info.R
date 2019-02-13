# CHNOSZ/protein.info.R

# calculate formulas and summarize properties of proteins
# pinfo: find rownumber in thermo$protein
# protein.length: lengths of the indicated proteins
# protein.formula: chemical makeup of the indicated proteins
# protein.obigt: perform group additivity calculations
# protein.basis: coefficients of basis species in formation reactions of [ionized] proteins [residues]
# protein.equil: step-by-step example of protein equilibrium calculation

pinfo <- function(protein, organism=NULL, residue=FALSE, regexp=FALSE) {
  # return the `protein` (possibly per residue) for:
  #   dataframe `protein`
  # return the rownumber(s) of thermo$protein for:
  #   character `protein`, e.g. LYSC_CHICK
  #   character `protein` and `organism`, e.g. 'LYSC', 'CHICK'
  # return the row(s) of thermo$protein (possibly per residue) for:
  #   numeric `protein` (the rownumber itself)
  t_p <- get("thermo", CHNOSZ)$protein
  if(is.data.frame(protein)) out <- protein
  if(is.numeric(protein)) {
    # drop NA matches to thermo$protein
    iproteins <- 1:nrow(t_p)
    protein[!protein %in% iproteins] <- NA
    # get amino acid counts
    out <- t_p[protein, ]
  }
  if(is.data.frame(protein) | is.numeric(protein)) {
    # compute per-residue counts if requested
    if(residue) out[, 5:25] <- out[, 5:25]/rowSums(out[, 6:25])
  } else {
    # search for protein by regular expression
    if(regexp) {
      iprotein <- grepl(protein, t_p$protein)
      iorganism <- iprotein
      if(!is.null(organism)) iorganism <- grepl(organism, t_p$organism)
      iprotein <- which(iprotein & iorganism)
    } else {
      # search for protein or protein_organism in thermo$protein
      t_p_names <- paste(t_p$protein, t_p$organism, sep="_")
      if(is.null(organism)) my_names <- protein
      else my_names <- paste(protein, organism, sep="_")
      iprotein <- match(my_names, t_p_names)
    }
    out <- iprotein
  }
  out
}

protein.formula <- function(protein, organism=NULL, residue=FALSE) {
  # return a matrix with chemical formulas of proteins
  aa <- pinfo(pinfo(protein, organism))
  rf <- group.formulas()
  out <- as.matrix(aa[, 5:25]) %*% as.matrix(rf)
  if(residue) out <- out / rowSums(aa[, 6:25])
  row.names(out) <- make.unique(paste(aa$protein, aa$organism, sep="_"))
  return(out)
}

protein.length <- function(protein, organism=NULL) {
  # calculate the length(s) of proteins
  aa <- pinfo(pinfo(protein, organism))
  # use rowSums on the columns containing amino acid counts
  pl <- as.numeric(rowSums(aa[, 6:25]))
  return(pl)
}

protein.obigt <- function(protein, organism=NULL, state=thermo()$opt$state) {
  # display and return the properties of
  # proteins calculated from amino acid composition
  aa <- pinfo(pinfo(protein, organism))
  # the names of the protein backbone groups depend on the state
  # [UPBB] for aq or [PBB] for cr
  if(state=="aq") bbgroup <- "UPBB" else bbgroup <- "PBB"
  # names of the AABB, sidechain and protein backbone groups
  groups <- c("AABB", colnames(aa)[6:25], bbgroup)
  # put brackets around the group names
  groups <- paste("[", groups, "]", sep="")
  # the rownumbers of the groups in thermo$obigt
  groups_state <- paste(groups, state)
  obigt <- get("thermo", CHNOSZ)$obigt
  obigt_state <- paste(obigt$name, obigt$state)
  igroup <- match(groups_state, obigt_state)
  # the properties are in columns 8-20 of thermo$obigt
  groupprops <- obigt[igroup, 8:20]
  # the elements in each of the groups
  groupelements <- i2A(igroup)
  # a function to work on a single row of aa
  eosfun <- function(aa) {
    # numbers of groups: chains [=AABB], sidechains, protein backbone
    nchains <- as.numeric(aa[, 5])
    length <- sum(as.numeric(aa[, 6:25]))
    npbb <- length - nchains
    ngroups <- c(nchains, as.numeric(aa[, 6:25]), npbb)
    # the actual adding and multiplying of thermodynamic properties
    # hmm. seems like we have to split up the multiplication/transposition
    # operations to get the result into multiple columns. 20071213
    eos <- t(data.frame(colSums(groupprops * ngroups)))
    # to get the formula, add up and round the group compositions 20090331
    f.in <- round(colSums(groupelements * ngroups), 3)
    # take out any elements that don't appear (sometimes S)
    f.in <- f.in[f.in!=0]
    # turn it into a formula
    f <- as.chemical.formula(f.in)
    # now the species name
    name <- paste(aa$protein, aa$organism, sep="_")
    # make some noise for the user
    message("protein.obigt: found ", appendLF=FALSE)
    message(name, " (", f, ", ", appendLF=FALSE)
    message(round(length, 3), " residues)")
    ref <- aa$ref
    header <- data.frame(name=name, abbrv=NA, formula=f, state=state, ref1=ref, ref2=NA, date=NA, stringsAsFactors=FALSE)
    eosout <- cbind(header, eos)
    return(eosout)
  }
  # loop over each row of aa
  out <- lapply(1:nrow(aa), function(i) eosfun(aa[i, ]))
  out <- do.call(rbind, out)
  rownames(out) <- NULL
  return(out)
}

protein.basis <- function(protein, T=25, normalize=FALSE) {
  # 20090902 calculate the coefficients of basis species in reactions
  # to form proteins (possibly per normalized by length) listed in protein
  # 20120528 renamed protein.basis from residue.info ...
  # what are the elemental compositions of the proteins
  aa <- pinfo(pinfo(protein))
  pf <- protein.formula(aa)
  # what are the coefficients of the basis species in the formation reactions
  sb <- species.basis(pf)
  # calculate ionization states if H+ is a basis species
  thermo <- get("thermo", CHNOSZ)
  iHplus <- match("H+", rownames(thermo$basis))
  if(!is.na(iHplus)) {
    pH <- -thermo$basis$logact[iHplus]
    Z <- ionize.aa(aa, T=T, pH=pH)[1, ]
    sb[, iHplus] <- sb[, iHplus] + Z
  }
  # compute per length-normalized coefficients if requested
  if(normalize) {
    # get lengths of proteins
    plen <- protein.length(aa)
    sb <- sb/plen
  }
  # return the result
  return(sb)
}

protein.equil <- function(protein, T=25, loga.protein=0, digits=4) {
  out <- character()
  mymessage <- function(...) {
    message(...)
    text <- paste(list(...), collapse = " ")
    out <<- c(out, text)
  }
  # show the individual steps in calculating metastable equilibrium among proteins
  mymessage("protein.equil: temperature from argument is ", T, " degrees C")
  TK <- convert(T, "K")
  # get the amino acid compositions of the proteins
  aa <- pinfo(pinfo(protein))
  # get some general information about the proteins
  pname <- paste(aa$protein, aa$organism, sep="_")
  plength <- protein.length(aa)
  # use thermo$basis to decide whether to ionize the proteins
  thermo <- get("thermo", CHNOSZ)
  ionize.it <- FALSE
  iword <- "nonionized"
  bmat <- basis.elements()
  if("H+" %in% rownames(bmat)) {
    ionize.it <- TRUE
    iword <- "ionized"
    pH <- -thermo$basis$logact[match("H+", rownames(bmat))]
    mymessage("protein.equil: pH from thermo$basis is ", pH)
  }
  # tell the user whose [Met] is in thermo$obigt
  info.Met <- info(info('[Met]', "aq"))
  mymessage("protein.equil: [Met] is from reference ", info.Met$ref1)
  ## first set of output: show results of calculations for a single protein
  mymessage("protein.equil [1]: first protein is ", pname[1], " with length ", plength[1])
  # standard Gibbs energies of basis species
  G0basis <- unlist(suppressMessages(subcrt(thermo$basis$ispecies, T=T, property="G")$out))
  # coefficients of basis species in formation reactions of proteins
  protbasis <- suppressMessages(protein.basis(aa, T=T))
  # sum of standard Gibbs energies of basis species in each reaction
  G0basissum <- colSums(t(protbasis) * G0basis)
  # standard Gibbs energies of nonionized proteins
  G0prot <- unlist(suppressMessages(subcrt(pname, T=T, property="G")$out))
  # standard Gibbs energy of formation reaction of nonionized protein, cal/mol
  G0protform <- G0prot - G0basissum
  mymessage("protein.equil [1]: reaction to form nonionized protein from basis species has G0(cal/mol) of ", signif(G0protform[1], digits))
  if(ionize.it) {
    # standard Gibbs energy of ionization of protein, cal/mol
    G0ionization <- suppressMessages(ionize.aa(aa, property="G", T=T, pH=pH))[1, ]
    mymessage("protein.equil [1]: ionization reaction of protein has G0(cal/mol) of ", signif(G0ionization[1], digits))
    # standard Gibbs energy of formation reaction of ionized protein, cal/mol
    G0protform <- G0protform + G0ionization
  }
  # standard Gibbs energy of formation reaction of non/ionized residue equivalents, dimensionless
  R <- 1.9872  # gas constant, cal K^-1 mol^-1
  G0res.RT <- G0protform/R/TK/plength
  mymessage("protein.equil [1]: per residue, reaction to form ", iword, " protein from basis species has G0/RT of ", signif(G0res.RT[1], digits))
  # coefficients of basis species in formation reactions of residues
  resbasis <- suppressMessages(protein.basis(aa, T=T, normalize=TRUE))
  # logQstar and Astar/RT
  logQstar <- colSums(t(resbasis) * - thermo$basis$logact)
  mymessage("protein.equil [1]: per residue, logQstar is ", signif(logQstar[1], digits))
  Astar.RT <- -G0res.RT - log(10)*logQstar
  mymessage("protein.equil [1]: per residue, Astar/RT = -G0/RT - 2.303logQstar is ", signif(Astar.RT[1], digits))
  if(!is.numeric(protein)) mymessage("protein.equil [1]: not comparing calculations with affinity() because 'protein' is not numeric")
  else {
    # for **Astar** we have to set the activities of the proteins to zero, not loga.protein!
    a <- suppressMessages(affinity(iprotein=protein, T=T, loga.protein=0))
    aAstar.RT <- log(10) * as.numeric(a$values) / plength
    mymessage("check it!       per residue, Astar/RT calculated using affinity() is ", signif(aAstar.RT[1], digits))
    if(!isTRUE(all.equal(Astar.RT, aAstar.RT, check.attributes=FALSE)))
      stop("Bug alert! The same value for Astar/RT cannot be calculated manually as by using affinity()")
  }
  if(length(pname)==1) mymessage("protein.equil [all]: all done... give me more than one protein for equilibrium calculations")
  else {
    ## next set of output: equilibrium calculations
    mymessage("protein.equil [all]: lengths of all proteins are ", paste(plength, collapse=" "))
    mymessage("protein.equil [all]: Astar/RT of all residue equivalents are ", paste(signif(Astar.RT, digits), collapse=" "))
    expAstar.RT <- exp(Astar.RT)
    sumexpAstar.RT <- sum(expAstar.RT)
    mymessage("protein.equil [all]: sum of exp(Astar/RT) of all residue equivalents is ", signif(sumexpAstar.RT, digits))
    # boltzmann distribution
    alpha <- expAstar.RT / sumexpAstar.RT    
    mymessage("protein.equil [all]: equilibrium degrees of formation (alphas) of residue equivalents are ", paste(signif(alpha, digits), collapse=" "))
    # check with equilibrate()
    if(is.numeric(protein)) {
      loga.equil.protein <- unlist(suppressMessages(equilibrate(a, normalize=TRUE))$loga.equil)
      # here we do have to convert from logarithms of activities of proteins to degrees of formation of residue equivalents
      a.equil.residue <- plength*10^loga.equil.protein
      ealpha <- a.equil.residue/sum(a.equil.residue)
      mymessage("check it!     alphas of residue equivalents from equilibrate() are ", paste(signif(ealpha, digits), collapse=" "))
      if(!isTRUE(all.equal(alpha, ealpha, check.attributes=FALSE)))
        stop("Bug alert! The same value for alpha cannot be calculated manually as by using equilibrate()")
    }
    # total activity of residues
    loga.residue <- log10(sum(plength * 10^loga.protein))
    mymessage("protein.equil [all]: for activity of proteins equal to 10^", signif(loga.protein, digits), ", total activity of residues is 10^", signif(loga.residue, digits))
    # equilibrium activities of residues
    loga.residue.equil <- log10(alpha*10^loga.residue)
    mymessage("protein.equil [all]: log10 equilibrium activities of residue equivalents are ", paste(signif(loga.residue.equil, digits), collapse=" "))
    # equilibrium activities of proteins
    loga.protein.equil <- log10(10^loga.residue.equil/plength)
    mymessage("protein.equil [all]: log10 equilibrium activities of proteins are ", paste(signif(loga.protein.equil, digits), collapse=" "))
    # check with equilibrate()
    if(is.numeric(protein)) {
      eloga.protein.equil <- unlist(suppressMessages(equilibrate(a, loga.balance=loga.residue, normalize=TRUE))$loga.equil)
      mymessage("check it!    log10 eq'm activities of proteins from equilibrate() are ", paste(signif(eloga.protein.equil, digits), collapse=" "))
      if(!isTRUE(all.equal(loga.protein.equil, eloga.protein.equil, check.attributes=FALSE)))
        stop("Bug alert! The same value for log10 equilibrium activities of proteins cannot be calculated manually as by using equilibrate()")
    }
  }
  return(out)
}

