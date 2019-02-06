## cross-checking sources
# the reference sources
ref.source <- thermo$refs$key
# sources in the primary thermodynamic database
# we omit the [S92] in "HDNB78 [S92]" etc.
tdata <- get("thermo")$obigt
ps1 <- gsub("\ .*", "", tdata$ref1)
ps2 <- gsub("\ .*", "", tdata$ref2)
# sources in the optional datafiles
tdata <- read.csv(system.file("extdata/OBIGT/DEW_aq.csv", package="CHNOSZ"), as.is=TRUE)
os1 <- gsub("\ .*", "", tdata$ref1)
os2 <- gsub("\ .*", "", tdata$ref2)
tdata <- read.csv(system.file("extdata/OBIGT/SLOP98.csv", package="CHNOSZ"), as.is=TRUE)
os3 <- gsub("\ .*", "", tdata$ref1)
os4 <- gsub("\ .*", "", tdata$ref2)
tdata <- read.csv(system.file("extdata/OBIGT/SUPCRT92.csv", package="CHNOSZ"), as.is=TRUE)
os5 <- gsub("\ .*", "", tdata$ref1)
os6 <- gsub("\ .*", "", tdata$ref2)
tdata <- read.csv(system.file("extdata/OBIGT/OldAA.csv", package="CHNOSZ"), as.is=TRUE)
os7 <- gsub("\ .*", "", tdata$ref1)
os8 <- gsub("\ .*", "", tdata$ref2)
# all of the thermodynamic data sources - some of them might be NA
obigt.source <- unique(c(ps1, ps2, os1, os2, os3, os4, os5, os6, os7, os8))
obigt.source <- obigt.source[!is.na(obigt.source)]
# these all produce character(0) if the sources are all accounted for
print("missing these sources for thermodynamic properties:")
print(unique(obigt.source[!(obigt.source %in% ref.source)]))
# determine if all the reference sources are cited
# this should produce character(0)
print("these sources are present but not cited:")
print(ref.source[!(ref.source %in% obigt.source)])
