# read.supcrt.R
# read thermodynamic data in SUPCRT format

# import values from supcrt files (e.g., slop98.dat)
# first version 20051105 jmd
read.supcrt <- function(file) {
  # the file must not have any comments in the data blocks
  # (the reference and data must be on a line by themselves)
  # 999999 will be converted to NA

  # to work correctly, a separator line i.e.
  #
  #  ******
  #
  # may have to be added after the aqueous species block,
  # before the table of number of species at the end

  # read the entire file
  tab <- scan(file,what='character')

  # function to identify separators
  # (long strings that begin with *)
  is.separator <- function(s) {
    if( nchar(s) > 5 & identical(substr(s,1,1),'*') ) return(TRUE)
    return(FALSE)
  }

  # function to find the first data position
  # for our desired state
  ifirst <- function(tab,state,istart=1,transitions=NULL) {
    for(i in istart:length(tab)) {
      if(is.separator(tab[i])) {
        if( identical(state,'aq') & identical(tab[i-1],'species') ) return(i+1)
        if( identical(state,'gas') & identical(tab[i-1],'gases' ) ) return(i+1)
        if( identical(state,'cr') & identical(transitions,0) 
          & identical(tab[i-3],'undergo') ) return(i+1)
        if( identical(state,'cr') & identical(transitions,1) 
          & identical(tab[i-3],'one') ) return(i+1)
        if( identical(state,'cr') & identical(transitions,2) 
          & identical(tab[i-3],'two') ) return(i+1)
        if( identical(state,'cr') & identical(transitions,3) 
          & identical(tab[i-3],'three') ) return(i+1)
      }
    }
  }

  name.lower <- function(name,abbrv,formula,state) {
    name <- as.character(name)
    abbrv <- as.character(abbrv)
    formula <- as.character(formula)
    for(i in 1:length(name)) {
      if( state=='aq' & (identical(name[i],formula[i]) 
        | identical(name[i],abbrv[i]) | identical(abbrv[i],formula[i]) ) ) next
      if( state=='cr' & (identical(name[i],formula[i]) 
        | identical(name[i],abbrv[i])) ) next
      name[i] <- tolower(name[i])
    }
    return(name)
  }
  
  strip <- function(s) {
    for(i in 1:length(s)) {
      si <- s[i]
      if(nchar(si)>3) {
        ss <- substr(si,nchar(si)-2,nchar(si))
        if(ss %in% c('(g)',',aq',',AQ','(s)','(-)','(+)','(0)','(1)','(2)') ) {
          s[i] <- substr(si,1,nchar(si)-3)
          if(identical(ss,'(-)')) s[i] <- paste(s[i],'-',sep='')
          if(identical(ss,'(+)')) s[i] <- paste(s[i],'+',sep='')
          if(identical(ss,'(1)')) s[i] <- paste(s[i],'+1',sep='')
          if(identical(ss,'(2)')) s[i] <- paste(s[i],'+2',sep='')
        }
      }
      if(nchar(si)>4) {
        ss <- substr(si,nchar(si)-3,nchar(si)) 
        if(ss %in% c('(+0)','(+1)','(+2)','(+3)','(+4)',
          '(-1)','(-2)','(-3)','(aq)','.dat') ) {
          s[i] <- substr(si,1,nchar(si)-4)
          if(identical(ss,'(-3)')) s[i] <- paste(s[i],'-3',sep='')
          if(identical(ss,'(-2)')) s[i] <- paste(s[i],'-2',sep='')
          if(identical(ss,'(-1)')) s[i] <- paste(s[i],'-1',sep='')
          if(identical(ss,'(+1)')) s[i] <- paste(s[i],'+1',sep='')
          if(identical(ss,'(+2)')) s[i] <- paste(s[i],'+2',sep='')
          if(identical(ss,'(+3)')) s[i] <- paste(s[i],'+3',sep='')
          if(identical(ss,'(+4)')) s[i] <- paste(s[i],'+4',sep='')
        }
      }
    }
    return(s)
  }

  source <- toupper(strip(file))
  
  read.cr.transitions <- function(tab,istart,transitions) {
    if(identical(transitions,0)) ni <- 14
    if(identical(transitions,1)) ni <- 21
    if(identical(transitions,2)) ni <- 28
    if(identical(transitions,3)) ni <- 35
    cr <- list()
    cr[[ni]] <- numeric()
    for(i in seq(istart,length(tab),ni)) {
      if(is.separator(tab[i])) break
      for(j in 1:ni) {
        cr[[j]] <- c(cr[[j]],tab[i+j-1])
      }
    }
    ghs <- data.frame( name=name.lower(strip(cr[[1]]),strip(cr[[3]]),
      strip(cr[[2]]),'cr'), abbrv=cr[[3]], formula=strip(cr[[2]]),
      state=rep('cr',length(cr[[1]])), source=strip(cr[[5]]),
      date=strip(cr[[6]]), Gf=cr[[7]], Hf=cr[[8]], S=cr[[9]] , stringsAsFactors=FALSE)
    eos <- data.frame( name=ghs$name, source=strip(cr[[5]]),
      a=cr[[11]], b=cr[[12]], c=cr[[13]], V=cr[[10]], T=cr[[14]] , stringsAsFactors=FALSE)
    trans <- data.frame(stringsAsFactors=FALSE)
    for(j in 0:transitions) {
      if(identical(transitions,0)) break
      it <- 1:7 + 7 + 7 * j
      t <- data.frame( name=ghs$name, Htr=cr[[it[1]]], Vtr=cr[[it[2]]], dP.dT=cr[[it[3]]],
        a=cr[[it[4]]], b=cr[[it[5]]], c=cr[[it[6]]], T=cr[[it[7]]] )
      trans <- rbind(trans,t)
    }
    return(list(ghs=ghs,eos=eos,transitions=trans,ilast=i))
  }

  read.cr <- function(tab) {
    icr0 <- ifirst(tab,'cr',1,0)
    cr0 <- read.cr.transitions(tab,icr0,0)
    cat(paste('read.supcrt:',nrow(cr0$ghs),
      'minerals that do not undergo phase transition\n'))
    icr1 <- ifirst(tab,'cr',cr0$ilast,1)
    cr1 <- read.cr.transitions(tab,icr1,1)
    cat(paste('read.supcrt:',nrow(cr1$ghs),
      'minerals that undergo one phase transition\n'))
    icr2 <- ifirst(tab,'cr',cr1$ilast,2)
    cr2 <- read.cr.transitions(tab,icr2,2)
    cat(paste('read.supcrt:',nrow(cr2$ghs),
      'minerals that undergo two phase transitions\n'))
    icr3 <- ifirst(tab,'cr',cr2$ilast,3)
    cr3 <- read.cr.transitions(tab,icr3,3)
    cat(paste('read.supcrt:',nrow(cr3$ghs),
      'minerals that undergo three phase transitions\n'))
    ghs <- rbind(cr0$ghs,cr1$ghs,cr2$ghs,cr3$ghs)
    eos <- rbind(cr0$eos,cr1$eos,cr2$eos,cr3$eos)
    transitions <- rbind(cr0$transition,cr1$transitions,cr2$transitions,cr3$transitions)
    return(list(ghs=ghs,eos=eos,transitions=transitions,ilast=cr3$ilast))
  }

  read.gas <- function(tab,istart) {
    gas <- list()
    gas[[14]] <- numeric()
    # populating the list
    for(i in seq(istart,length(tab),14)) {
      if(is.separator(tab[i])) break
      for(j in 1:14) gas[[j]] <- c(gas[[j]],tab[i+j-1])
    }
    # processing the list
    ghs <- data.frame( name=tolower(gas[[2]]), 
      abbrv=strip(gas[[3]]), formula=strip(gas[[3]]),
      state=rep('gas',length(gas[[1]])), source=strip(gas[[5]]),
      date=strip(gas[[6]]),Gf=gas[[7]], Hf=gas[[8]], S=gas[[9]] )
    eos <- data.frame( name=ghs$name, source=strip(gas[[5]]),
      a=gas[[11]], b=gas[[12]], c=gas[[13]], V=gas[[10]], T=gas[[14]] )
    return(list(ghs=ghs,eos=eos,ilast=i))
  }

  read.aq <- function(tab,istart) {
    aq <- list()
    aq[[17]] <- numeric()
    for(i in seq(istart,length(tab),17)) {
      if(is.separator(tab[i])) break
      for(j in 1:17) aq[[j]] <- c(aq[[j]],tab[i+j-1])
    }
    ghs <- data.frame( name=name.lower(strip(aq[[1]]),strip(aq[[3]]),
      strip(aq[[2]]),'aq'), abbrv=strip(aq[[3]]), formula=strip(aq[[2]]),
      state=rep('aq',length(aq[[1]])), source=strip(aq[[5]]),
      date=strip(aq[[6]]), Gf=aq[[7]], Hf=aq[[8]], S=aq[[9]] )
    eos <- data.frame( name=ghs$name, source=strip(aq[[5]]),
      a1=aq[[10]], a2=aq[[11]], a3=aq[[12]], a4=aq[[13]], 
      c1=aq[[14]], c2=aq[[15]], omega=aq[[16]] )
    return(list(ghs=ghs,eos=eos,ilast=i))
  }

  cr <- read.cr(tab)
  igas <- ifirst(tab,'gas',cr$ilast)
  gas <- read.gas(tab,igas)
  cat(paste('read.supcrt:',nrow(gas$ghs),'gases\n'))
  iaq <- ifirst(tab,'aq',gas$ilast)
  aq <- read.aq(tab,iaq)
  cat(paste('read.supcrt:',nrow(aq$ghs),'aqueous species\n'))

  supcrt <- list(ghs=rbind(cr$ghs,gas$ghs,aq$ghs),
    eos.cr=cr$eos,eos.gas=gas$eos,eos.aq=aq$eos,transitions=cr$transitions)
  return(supcrt)

  # some numerical conversions
  as.speciate <- function(supcrt) {
    for(i in 1:ncol(supcrt)) {
      supcrt[,i] <- as.numeric(as.character(supcrt[,i]))
      supcrt[supcrt[,i]==999999,i] <- NA
      cname <- colnames(supcrt)[i]
      if(cname=='b') supcrt[,i] <- supcrt[,i] * 10^(-3)
      if(cname=='c') supcrt[,i] <- supcrt[,i] * 10^(5)
      if(cname=='a1') supcrt[,i] <- supcrt[,i] * 10^(-1)
      if(cname=='a2') supcrt[,i] <- supcrt[,i] * 10^(2)
      if(cname %in% c('a4','c2')) supcrt[,i] <- supcrt[,i] * 10^(4)
      if(cname=='omega') supcrt[,i] <- supcrt[,i] * 10^(5)
    }
    return(supcrt)
  }

  supcrt$ghs[,6:8] <- as.speciate(supcrt$ghs[,6:8])
  supcrt$eos.cr[,3:7] <- as.speciate(supcrt$eos.cr[,3:7])
  supcrt$transitions[,2:8] <- as.speciate(supcrt$transitions[,2:8])
  supcrt$eos.gas[,3:7] <- as.speciate(supcrt$eos.gas[,3:7])
  supcrt$eos.aq[,3:9] <- as.speciate(supcrt$eos.aq[,3:9])

  return(supcrt)
}

