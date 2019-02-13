\encoding{UTF-8}
\name{buffer}
\alias{buffer}
\alias{mod.buffer}
\title{Calculating Buffered Chemical Activities}
\description{
Calculate values of activity or fugacity of basis species buffered by an assemblage of one or more species.
}

\usage{
  mod.buffer(name, species = NULL, state = thermo()$opt$state, 
    logact = -3)
}

\arguments{
  \item{name}{character, name of buffer to add to or find in \code{thermo$buffers}.}
  \item{species}{character, names or formulas of species in a buffer.}
  \item{state}{character, physical states of species in buffer.}
  \item{logact}{numeric, logarithms of activities of species in buffer.}
}

\details{
A buffer is treated here as assemblage of one or more species whose presence constrains values of the chemical activity (or fugacity) of one or more basis species.
To perform calculations for buffers use \code{\link{basis}} to associate the name of the buffer with one or more basis species.
After this, calls to \code{affinity} will invoke the required calculations.
The calculated values of the buffered activites can be retrieved by setting \code{return.buffer} to \code{TRUE} (in \code{\link{affinity}}).
The maximum number of buffered chemical activities possible for any buffer is equal to the number of species in the buffer; however, the user may then elect to work with the values for only one or some of the basis species calculated with the buffer. 

The identification of a conserved basis species (or other reaction balancing rule) is required in calculations for buffers of more than one species.
For example, in the pyrite-pyrrhotite-magnetite buffer (\eqn{\mathrm{FeS_2}}{FeS2}-\eqn{\mathrm{FeS}}{FeS}-\eqn{\mathrm{Fe_3O_4}}{Fe3O4}) a basis species common to each species is one representing \eqn{Fe}{Fe}.
Therefore, when writing reactions between the species in this buffer \eqn{Fe}{Fe} is conserved while \eqn{\mathrm{H_2S}}{H2S} and \eqn{\mathrm{O2}}{O2} are the variables of interest.
The calculation for buffers attempts to determine which of the available basis species qualifies as a conserved quantity.
This can be overriden with \code{balance}.
The default value of \code{balance} is \samp{PBB}, which instructs the function to use the protein backbone group as the conserved quantity in buffers consisting of proteins, but has no overriding effect on the computations for buffers without proteins.

To view the available buffers, print the \code{\link{thermo}$buffer} object.
Buffer definitions can be added to this dataframe with \code{mod.buffer}.
It is possible to set the logarithms of activities of the species in the buffer through the \code{logact} argument; if this is missing unit activity is assigned to crystalline species in buffer, otherwise (for aqueous species) the default value of activity is 10\eqn{^{-3}}{^-3}.
If \code{name} identifies an already defined buffer, this function modifies the logarithms of activities or states of species in that buffer, optionally restricted to only those species given in \code{species}.

It is possible to assign different buffers to different basis species, in which case the order of their calculation depends on their order in \code{thermo$buffers}.
This function is compatible with systems of proteins, but note that for buffers \emph{made} of proteins the buffer calculations presently use whole protein formulas (instead of residue equivalents) and consider nonionized proteins only.
}

\seealso{
\code{\link{protein}} for an example using a buffer made of proteins.
}

\examples{
\dontshow{reset()}
## list the buffers
thermo()$buffers
# another way to do it, for a specific buffer
print(mod.buffer("PPM"))

## buffer made of one species
# calculate the activity of CO2 in equilibrium with
# (a buffer made of) acetic acid at a given activity
basis("CHNOS")
basis("CO2","AC")
# what activity of acetic acid are we using?
print(mod.buffer("AC"))
# return the activity of CO2
(logaCO2 <- affinity(return.buffer=TRUE)$CO2)
stopifnot(all.equal(logaCO2, -7.05752136))
# as a function of oxygen fugacity
affinity(O2=c(-85,-70,4),return.buffer=TRUE)
# as a function of logfO2 and temperature
affinity(O2=c(-85,-70,4),T=c(25,100,4),return.buffer=TRUE)
# change the activity of species in the buffer
mod.buffer("AC",logact=-10)
affinity(O2=c(-85,-70,4),T=c(25,100,4),return.buffer=TRUE)
# see below for a different strategy using the
# 'type' argument of diagram

## buffer made of three species
## Pyrite-Pyrrhotite-Magnetite (PPM)
# specify basis species and initial activities
basis(c("FeS2","H2S","O2","H2O"),c(0,-10,-50,0))
# note that the affinity of formation of pyrite,
# which corresponds to FeS2 in the basis, is zero
species(c("pyrite","pyrrhotite","magnetite"))
affinity(T=c(200,400,11),P=2000)$values
# setup H2S and O2 to be buffered by PPM
basis(c("H2S","O2"),c("PPM","PPM"))
# inspect values of H2S activity and O2 fugacity
affinity(T=c(200, 400, 11), P=2000, return.buffer=TRUE, exceed.Ttr=TRUE)
# now, the affinities of formation reactions of 
# species in the buffer are all equal to zero
print(a <- affinity(T=c(200, 400, 11), P=2000, 
  exceed.Ttr=TRUE)$values)
for(i in 1:length(a)) stopifnot(isTRUE(
  all.equal(as.numeric(a[[i]]),rep(0,length(a[[i]])))))

## buffer made of one species: show values of logfO2 on an 
## Eh-pH diagram; after Garrels, 1960, Figure 6
basis("CHNOSe")
# here we will buffer the activity of the electron by O2
mod.buffer("O2","O2","gas",999)
basis("e-","O2")
# start our plot, then loop over values of logfO2
thermo.plot.new(xlim=c(0,14),ylim=c(-0.8,1.2),
  xlab="pH",ylab=axis.label("Eh"))
# the upper and lower lines correspond to the upper
# and lower stability limits of water
logfO2 <- c(0,-20,-40,-60,-83.1)
for(i in 1:5) {
  # update the logarithm of fugacity (logact) of O2 in the buffer
  mod.buffer("O2","O2","gas",logfO2[i])
  # get the values of the logarithm of activity of the electron
  a <- affinity(pH=c(0,14,15),return.buffer=TRUE)
  # convert values of pe (-logact of the electron) to Eh
  Eh <- convert(-as.numeric(a$`e-`),"Eh")
  lines(seq(0,14,length.out=15),Eh)
  # add some labels
  text(seq(0,14,length.out=15)[i*2+2],Eh[i*2+2],
    paste("logfO2=",logfO2[i],sep=""))
}
title(main=paste("Relation between logfO2(g), Eh and pH at\n",
  "25 degC and 1 bar. After Garrels, 1960"))

## buffer made of two species
# conditions for metastable equilibrium among 
# CO2 and acetic acid. note their starting activities:
print(mod.buffer("CO2-AC")) 
basis("CHNOS")
basis("O2","CO2-AC")
affinity(return.buffer=TRUE)  # logfO2 = -75.94248
basis("CO2",123)  # what the buffer reactions are balanced on
affinity(return.buffer=TRUE)  # unchanged
# consider more oxidizing conditions
mod.buffer("CO2-AC",logact=c(0,-10))
affinity(return.buffer=TRUE)

# one can solve for the logarithm of activity of a 
# basis species using the 'type' argument of diagram
basis("CHNOS")
basis("CO2", 999)
species("acetic acid", -3)
a <- affinity(O2=c(-85, -70, 4), T=c(25, 100, 4))
# write a title with formulas and subscripts
lCO2 <- axis.label("CO2")
main <- substitute(a~~b~~c,list(a=lCO2, b="buffered by",
  c="acetic acid"))
d <- diagram(a, type="CO2", main=main)
species(1, -10)
a <- affinity(O2=c(-85, -70, 4), T=c(25, 100, 4))
d <- diagram(a, type="CO2", add=TRUE, lty=2)
# add a legend
lAC <- expr.species("CH3COOH", log=TRUE)
ltext <- c(as.expression(lAC), -3, -10)
lty <- c(NA, 1, 2)
legend("topright", legend=ltext, lty=lty, bg="white")
# do return.buffer and diagram(type=...) give the same results?
and <- as.numeric(d$plotvals[[1]])
basis("CO2", "AC")
mod.buffer("AC", logact=-10)
a.buffer <- affinity(O2=c(-85, -70, 4), T=c(25, 100, 4),
  return.buffer=TRUE)
ana <- as.numeric(unlist(a.buffer[[1]]))
stopifnot(all.equal(ana, and))
}

\references{
Garrels, R. M. (1960) \emph{Mineral Equilibria}. Harper & Brothers, New York, 254 p. \url{http://www.worldcat.org/oclc/552690}
}

\concept{Extended workflow}
