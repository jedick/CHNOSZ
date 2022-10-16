# CHNOSZ/util.fasta.R
# read and manipulate FASTA sequence files

read.fasta <- function(file, iseq=NULL, ret="count", lines=NULL, ihead=NULL,
  start=NULL, stop=NULL, type="protein", id=NULL) {
  # read sequences from a fasta file
  # some of the following code was adapted from 
  # read.fasta in package seqinR
  # value of 'iseq' is what sequences to read (default is all)
  # value of 'ret' determines format of return value:
  #   count: amino acid composition (same columns as thermo()$protein, can be used by add.protein)
  #        or nucleic acid base composition (A-C-G-T)
  #   seq: amino acid sequence
  #   fas: fasta entry
  # value of 'id' is used for 'protein' in output table,
  #   otherwise ID is parsed from FASTA header (can take a while)
  
  # check if the file is in an archive (https://github.com/jimhester/archive)
  if(inherits(file, "archive_read")) {
    is.archive <- TRUE
    filebase <- gsub("]", "", basename(summary(file)$description))
  } else {
    is.archive <- FALSE
    filebase <- basename(file)
  }
  if(is.null(lines)) {
    message("read.fasta: reading ", filebase, " ... ", appendLF=FALSE)
    is.nix <- Sys.info()[[1]]=="Linux"
    if(is.archive) {
      # we can't use scan here?
      lines <- readLines(file)
    } else if(is.nix) {
      # retrieve contents using system command (seems slightly faster even than scan())
      # figure out whether to use 'cat', 'zcat' or 'xzcat'
      suffix <- substr(file,nchar(file)-2,nchar(file))
      if(suffix==".gz") mycat <- "zcat"
      else if(suffix==".xz") mycat <- "xzcat"
      else mycat <- "cat"
      lines <- system(paste(mycat,' "',file,'"',sep=""),intern=TRUE)
    } else lines <- scan(file, what=character(), sep="\n", quiet=TRUE)
  }
  nlines <- length(lines)
  message(nlines, " lines ... ", appendLF=FALSE)
  if(is.null(ihead)) ihead <- which(substr(lines,1,1)==">")
  message(length(ihead), " sequences")
  linefun <- function(i1,i2) lines[i1:i2]
  # identify the lines that begin and end each sequence
  begin <- ihead + 1
  end <- ihead - 1
  end <- c(end[-1], nlines)
  # use all or selected sequences
  if(is.null(iseq)) iseq <- seq_along(begin)
  # just return the lines from the file
  if(ret=="fas") {
    iline <- numeric()
    for(i in iseq) iline <- c(iline,(begin[i]-1):end[i])
    return(lines[iline])
  }
  # get each sequence from the begin to end lines
  seqfun <- function(i) paste(linefun(begin[i],end[i]),collapse="")
  sequences <- lapply(iseq, seqfun)
  # organism name is from file name
  # (basename minus extension)
  bnf <- strsplit(filebase,split=".",fixed=TRUE)[[1]][1]
  organism <- bnf
  # protein/gene name is from header line for entry
  # (strip the ">" and go to the first space)
  missid <- missing(id)
  if(is.null(id)) id <- as.character(lapply(iseq, function(j) {
    # get the text of the line
    f1 <- linefun(ihead[j],ihead[j])
    # stop if the first character is not ">"
    # or the first two charaters are "> "
    if(substr(f1,1,1)!=">" | length(grep("^> ",f1)>0))
      stop(paste("file",filebase,"line",j,"doesn't begin with FASTA header '>'."))
    # discard the leading '>'
    f2 <- substr(f1, 2, nchar(f1))
    # keep everything before the first space
    return(strsplit(f2," ")[[1]][1])
  } ))
  if(ret=="count") {
    counts <- count.aa(sequences, start, stop, type)
    ref <- abbrv <- NA
    chains <- 1
    if(type=="protein") {
      colnames(counts) <- aminoacids(3)
      # 20090507 made stringsAsFactors FALSE
      out <- cbind(data.frame(protein=id, organism=organism,
        ref=ref, abbrv=abbrv, chains=chains, stringsAsFactors=FALSE), counts)
      # 20170117 extra processing for files from UniProt
      isUniProt <- grepl("\\|......\\|.*_", out$protein[1])
      if(isUniProt & missid) {
        p1 <- sapply(strsplit(out$protein, "\\|"), "[", 1)
        p2 <- sapply(strsplit(out$protein, "\\|"), "[", 2)
        p3 <- sapply(strsplit(out$protein, "\\|"), "[", 3)
        out$abbrv <- sapply(strsplit(p3, "_"), "[", 1)
        out$organism <- sapply(strsplit(p3, "_"), "[", 2)
        out$protein <- paste0(p1, "|", p2)
      }
      out
    } else if(type %in% c("DNA", "RNA")) {
      cbind(data.frame(gene=id, organism=organism,
        ref=ref, abbrv=abbrv, chains=chains, stringsAsFactors=FALSE), counts)
    }
  } else return(sequences)
}

uniprot.aa <- function(protein, start=NULL, stop=NULL) {
  # download protein sequence information from UniProt
  iprotein <- numeric()
  # construct the initial URL
  proteinURL <- paste("https://www.uniprot.org/uniprot/", protein, sep="")
  message("uniprot.aa: trying ", proteinURL, " ...", appendLF=FALSE)
  # try loading the URL, hiding any warnings
  oldopt <- options(warn=-1)
  URLstuff <- try(readLines(proteinURL),TRUE)
  options(oldopt)
  if(inherits(URLstuff, "try-error")) {
    message(" ::: FAILED :::")
    return(NA)
  }
  # 20091102: look for a link to a fasta file
  link <- grep("/uniprot/.*fasta", URLstuff)
  if(length(link) > 0) linkline <- URLstuff[[link[1]]]
  else {
    message(" ::: FAILED :::")
    return(NA)
  }
  # extract accession number from the link
  linkhead <- strsplit(linkline, ".fasta", fixed=TRUE)[[1]][1]
  accession.number <- tail(strsplit(linkhead, "/uniprot/", fixed=TRUE)[[1]], 1)
  message(" accession ", accession.number, " ...")
  # now download the fasta file
  fastaURL <- paste("https://www.uniprot.org/uniprot/", accession.number, ".fasta", sep="")
  URLstuff <- readLines(fastaURL)
  # get the header information / show  the user
  header <- URLstuff[[1]]
  header3 <- strsplit(header, "|", fixed=TRUE)[[1]][3]
  headerP_O <- strsplit(header3, " ")[[1]][1]
  header.id <- strsplit(header, headerP_O)[[1]][1]
  header.id <- substr(header.id, 2, nchar(header.id)-1)
  header.organism <- strsplit(headerP_O, "_")[[1]][2]
  message(paste0(header), appendLF=FALSE)
  # 20130206 use read.fasta with lines, start, stop arguments
  aa <- read.fasta(file="", lines=URLstuff, start=start, stop=stop)
  message(" (length ", sum(aa[1, 6:25]), ")", sep="")
  aa$protein <- header.id
  aa$organism <- header.organism
  return(aa)
}

count.aa <- function(seq, start=NULL, stop=NULL, type="protein") {
  # count amino acids or DNA bases in one or more sequences given as elements of the list seq
  if(type=="protein") letts <- aminoacids(1)
  else if(type=="DNA") letts <- c("A", "C", "G", "T")
  else if(type=="RNA") letts <- c("A", "C", "G", "U")
  else stop(paste("unknown sequence type", type))
  # the numerical positions of the letters in alphabetical order (i.e. for amino acids, same order as in thermo()$protein)
  ilett <- match(letts, LETTERS)
  # the letters A-Z represented by raw values
  rawAZ <- charToRaw("ABCDEFGHIJKLMNOPQRSTUVWXYZ")
  # to count the letters in each sequence
  countfun <- function(seq, start, stop) {
    # get a substring if one or both of start or stop are given
    # if only one of start or stop is given, get a default value for the other
    if(!is.null(start)) {
      if(is.null(stop)) stop <- nchar(seq)
      seq <- substr(seq, start, stop)
    } else if(!is.null(stop)) {
      seq <- substr(seq, 1, stop)
    }
    ## the actual counting ...
    #nnn <- table(strsplit(toupper(seq), "")[[1]])
    # ... replaced with C version 20180217
    counts <- .C(C_count_letters, seq, integer(26))[[2]]
    # which is equivalent to this R code:
    #rawseq <- charToRaw(toupper(seq))
    #counts <- sapply(rawAZ, function(x) sum(rawseq == x))
    return(counts)
  }
  # counts for each sequence
  counts <- lapply(seq, countfun, start, stop)
  counts <- do.call(rbind, counts)
  # check for letters that aren't in our alphabet
  ina <- colSums(counts[, -ilett, drop=FALSE]) > 0
  if(any(ina)) {
    message(paste("count.aa: unrecognized letter(s) in", type, "sequence:", paste(LETTERS[-ilett][ina], collapse=" ")))
  }
  counts <- counts[, ilett, drop=FALSE]
  # clean up row/column names
  colnames(counts) <- letts
  rownames(counts) <- 1:nrow(counts)
  return(counts)
}

