# CHNOSZ/buffer.R
# Calculate chemical activities of buffered species
# 20061102 jmd

mod.buffer <- function(name, species = NULL, state = "cr", logact = 0) {
  # 20071102 add or change a buffer system
  thermo <- get("thermo", CHNOSZ)
  if(is.null(species)) {
    iname <- which(name == thermo$buffer$name)
    if(length(iname) > 0) species <- thermo$buffer$species[iname]
    else species <- character()
  }
  ls <- length(species)
  if(ls < length(name) | ls < length(state) | ls < length(logact))
    stop('species must be at least as long as the other arguments')
  if(length(name) != ls) name <- rep(name, length.out = ls)
  add <- TRUE
  if(TRUE %in% (name %in% thermo$buffer$name)) {
    add <- FALSE
    imod <- which(thermo$buffer$name %in% name & thermo$buffer$species %in% species)
    if(length(imod)>0) {
      if(state[1] == '') {
        thermo$buffer <- thermo$buffer[-imod, ]
        assign("thermo", thermo, CHNOSZ)
        message(paste('mod.buffer: removed ', c2s(species), ' in ', c2s(unique(name)), ' buffer', sep = ''))
      } else {
        if(missing(state)) state <- thermo$buffer$state[imod]
        if(missing(logact)) logact <- thermo$buffer$logact[imod]
        if(length(state) != ls) state <- rep(state,length.out = ls)
        if(length(logact) != ls) logact <- rep(logact,length.out = ls)
        state.old <- thermo$buffer$state[imod]
        logact.old <- thermo$buffer$logact[imod]
        thermo$buffer$state[imod] <- state
        thermo$buffer$logact[imod] <- logact
        assign("thermo", thermo, CHNOSZ)
        if(identical(state.old, state) & identical(logact.old, logact)) {
          message(paste('mod.buffer: nothing changed for ', c2s(species), ' in ', c2s(unique(name)), ' buffer', sep = ''))
        } else {
          message(paste('mod.buffer: changed state and/or logact of ', c2s(species), ' in ', c2s(unique(name)), ' buffer', sep = ''))
        }
      }
    } else {
      add <- TRUE
    }
  } 
  if(add) {
    if(state[1] == '') state <- rep(thermo$opt$state, length.out = ls)
    t <- data.frame(name = name, species = species, state = state, logact = logact)
    thermo$buffer <- rbind(thermo$buffer, t)
    assign("thermo", thermo, CHNOSZ)
    message(paste('mod.buffer: added',c2s(unique(name))))
  }
  return(invisible(thermo$buffer[thermo$buffer$name %in% name, ]))
}

### Unexported functions ###

buffer <- function(logK = NULL, ibasis = NULL, logact.basis = NULL, is.buffer = NULL, balance = 'PBB') {
  thermo <- get("thermo", CHNOSZ)
  # If logK is NULL load the buffer species,
  # otherwise perform buffer calculations
  if(is.null(logK)) {
    # Load the buffer species
    buffers <- unique(as.character(thermo$basis$logact)[!can.be.numeric(as.character(thermo$basis$logact))])
    ispecies.new <- list()
    for(k in 1:length(buffers)) {
      ibasis <- which(thermo$basis$logact == buffers[k])
      ispecies <- numeric()
      for(i in 1:length(ibasis)) {
        ib <- as.character(thermo$buffer$name) == as.character(thermo$basis$logact[ibasis[i]])
        species <- as.character(thermo$buffer$species)[ib]
        state <- as.character(thermo$buffer$state)[ib]
        #ibuff <- info(species,state,quiet = TRUE)
        ispecies <- c(ispecies, species(species, state, index.return = TRUE, add = TRUE))
      }
      ispecies.new <- c(ispecies.new, list(ispecies))
      # Make sure to set the activities
      species(ispecies, thermo$buffer$logact[ib])
    }
    names(ispecies.new) <- buffers
    return(ispecies.new)
  }

  # Sometimes (e.g. PPM) the buffer species are identified multiple
  # times, causing problems for square matrices
  # Make the species appear singly
  is.buffer <- unique(is.buffer)
  bufbasis <- species.basis(thermo$species$ispecies[is.buffer])
  bufname <- thermo$basis$logact[ibasis[1]]
  basisnames <- rownames(thermo$basis)
  are.proteins <- grep('_', as.character(thermo$species$name[is.buffer]))
  if((length(are.proteins) > 0 & balance == 'PBB') | balance == 1) {
    if(balance == 1) {
      basisnames <- c('product', basisnames)
      nb <- rep(1, nrow(bufbasis))
      bufbasis <- cbind(data.frame(product = nb), bufbasis)
    } else {
      basisnames <- c('PBB', basisnames)
      # Prepend a PBB column to bufbasis and increment ibasis by 1
      nb <- as.numeric(protein.length(thermo$species$name[is.buffer]))
      bufbasis <- cbind(data.frame(PBB = nb), bufbasis)
    }
    ibasis <- ibasis + 1
    # Make logact.basis long enough
    tl <- length(logact.basis)
    logact.basis[[tl+1]] <- logact.basis[[tl]]
    # Rotate the entries so that the new one is first
    ilb <- c(tl+1, 1:tl)
    logact.basis <- logact.basis[ilb]
  }
  # Say hello
  #message(paste("buffer: '",bufname,"', of ",length(is.buffer),' species, ',length(ibasis),' activity(s) requested.',sep = ''))
  ibasisrequested <- ibasis
  # Check and maybe add to the number of buffered activities
  ibasisadded <- numeric()
  if( (length(ibasis)+1) != length(is.buffer) & length(is.buffer) > 1) {
    # Try to add buffered activities the user didn't specify
    # (e.g. H2S in PPM buffer if only H2 was requested)
    for(i in 1:(length(is.buffer) - (length(ibasis)+1))) {
      newbasis <- NULL
      # We want to avoid any basis species that might be used as the conservant
      # look for additional activities to buffer ... do columns in reverse 
      for(j in ncol(bufbasis):1) {
        if(j %in% ibasis) next
        if(FALSE %in% (bufbasis[,j] == 0)) {
          newbasis <- j
          break
        }
      }
      if(!is.null(newbasis)) {
        ibasis <- c(ibasis, newbasis)
        ibasisadded <- c(ibasisadded, newbasis)
      } else {
        stop('can not find enough buffered basis species for ', thermo$basis$logact[ibasis[1]], '.', sep = '')
      }
    }
  } 
  # And the leftovers
  #xx <- as.data.frame(bufbasis[,-ibasis])
  # The final buffered activity: the would-be conserved component
  newbasis <- NULL
  if(length(is.buffer) > 1) {
    # First try to get one that is present in all species
    for(i in ncol(bufbasis):1) {
      if(i %in% ibasis) next
      if(!TRUE %in% (bufbasis[, i] == 0)) newbasis <- i
    }
    # Or look for one that is present at all
    if(is.null(newbasis)) for(i in ncol(bufbasis):1) {
      if(i %in% ibasis) next
      if(FALSE %in% (bufbasis[, i] == 0)) newbasis <- i
    }
    if(!is.null(newbasis)) {
      ibasis <- c(ibasis,newbasis)
      #message(paste('buffer: the conserved activity is ', basisnames[newbasis], '.', sep = ''))
      #thermo$basis$logact[newbasis] <<- thermo$basis$logact[ibasis[1]]
    }
    else stop('no conserved activity found in your buffer (not enough basis species?)!')
  }
  if(is.null(newbasis)) context <- '' else context <- paste(', ', basisnames[newbasis], ' (conserved)', sep = '')
  reqtext <- paste(c2s(basisnames[ibasisrequested]), ' (active)', sep = '')
  if(length(ibasisadded) == 0) addtext <- '' else addtext <- paste(', ', c2s(basisnames[ibasisadded]), sep = '')
  message(paste("buffer: '", bufname, "' for activity of ", reqtext, addtext, context, sep = ""))
  # There could still be stuff here (over-defined system?)
  xx <- bufbasis[, -ibasis, drop = FALSE]
  # For the case when all activities are buffered
  if(ncol(xx) == 0) xx <- data.frame(xx = 0)
  # Our stoichiometric matrix - should be square
  A <- as.data.frame(bufbasis[,ibasis])
  # Determine conservation coefficients
  # Values for the known vector
  B <- list()
  for(i in 1:length(is.buffer)) {
    b <- -logK[[is.buffer[i]]] + thermo$species$logact[is.buffer[i]]
    if(ncol(xx) > 0) {
      if(is.list(xx)) xxx <- xx[[1]] else xxx <- xx
      if(ncol(xx) == 1 & identical(as.numeric(xxx), 0)) {
        # Do nothing
      } else {
        for(j in 1:ncol(xx)) {
          bs <- xx[i, j] * logact.basis[[match(colnames(xx)[j], basisnames)]]
          if(!is.matrix(bs)) bs <- matrix(bs, byrow = TRUE, nrow = nrow(as.data.frame(logact.basis[[1]])))
          if(ncol(bs) == 1) b <- matrix(b)
          b <- b - bs
        }
      }
    }
    # Force this to be matrix so indexing works at a single point (no variables defined in affinity()) 20201102
    B[[i]] <- as.matrix(b)
  }
  # A place to put the results
  X <- rep(B[1],length(ibasis))
  for(i in 1:nrow(B[[1]])) {
    for(j in 1:ncol(B[[1]])) {
      b <- numeric()
      for(k in 1:length(B)) b <- c(b, B[[k]][i, j])
      AAA <- A
      # solve for the activities in the buffer
      t <- solve(AAA, b)
      for(k in 1:length(ibasis))
        X[[k]][i, j] <- t[k]
    }
  }
  # Store results
  for(i in 1:length(ibasis)) {
    if(ncol(X[[i]]) == 1) X[[i]] <- as.numeric(X[[i]])
    else if(nrow(X[[i]]) == 1) X[[i]] <- as.matrix(X[[i]], nrow = 1)
    logact.basis[[ibasis[i]]] <- X[[i]]
  }
  names(logact.basis) <- basisnames
  
  return(list(ibasis = ibasis, logact.basis = logact.basis))   
}
