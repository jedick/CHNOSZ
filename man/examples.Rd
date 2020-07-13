\encoding{UTF-8}
\name{examples}
\alias{examples}
\alias{demos}
\title{Run Examples from the Documentation}
\description{
  Run the examples contained in each of the documentation topics. }

\arguments{
  \item{save.png}{logical, generate PNG image files for the plots?}
  \item{which}{character, which example to run}
}

\usage{
  examples(save.png = FALSE)
  demos(which = c("sources", "protein.equil", "affinity", "NaCl",
    "density", "ORP", "findit", "ionize", "buffer", "protbuff",
    "glycinate", "mosaic", "copper", "arsenic", "solubility", "gold",
    "contour", "sphalerite", "wjd", "bugstab", "Shh", "saturation",
    "adenine", "DEW", "lambda", "TCA", "aluminum", "AkDi",
    "comproportionation"),
    save.png=FALSE)
}

\details{

\code{examples} runs all the examples in the help pages for the package. 
\code{\link{example}} is called for each topic with \code{ask} set to \code{FALSE} (so all of the figures are shown without prompting the user).

\code{demos} runs all the \code{\link{demo}s} in the package.
The demo(s) to run is/are specified by \code{which}; the default is to run them in the order of the list below.

  \tabular{ll}{
    \code{sources} \tab Cross-check the reference list with the thermodynamic database \cr
    \code{protein.equil} \tab Chemical activities of two proteins in metastable equilibrium (Dick and Shock, 2011) \cr
    \code{affinity} \tab Affinities of metabolic reactions and amino acid synthesis (Amend and Shock, 1998, 2001) \cr
    \code{NaCl} \tab Equilibrium constant for aqueous NaCl dissociation (Shock et al., 1992) \cr
    \code{density} \tab Density of \H2O, inverted from IAPWS-95 equations (\code{\link{rho.IAPWS95}}) \cr
    \code{ORP} \tab Temperature dependence of oxidation-reduction potential for redox standards \cr
    \code{findit} \tab Minimize the standard deviation of logarithms of activities of sulfur species \cr
    \code{ionize} \tab ionize.aa(): contour plots of net charge and ionization properties of LYSC_CHICK \cr
    \code{buffer} \tab Minerals and aqueous species as buffers of hydrogen fugacity (Schulte and Shock, 1995) \cr
    \code{protbuff} \tab Chemical activities buffered by thiol peroxidases or sigma factors \cr
    \code{glycinate} \tab Metal-glycinate complexes (Shock and Koretsky, 1995; Azadi et al., 2019) \cr
    \code{mosaic} \tab Eh-pH diagram with two sets of changing basis species (Garrels and Christ, 1965) \cr
    \code{copper} \tab Another example of \code{\link{mosaic}}: complexation of Cu with glycine (Aksu and Doyle, 2001) \cr
    \code{arsenic} \tab Another example of \code{\link{mosaic}}: Eh-pH diagram for the system As-O-H-S (Lu and Zhu, 2011) \cr
    \code{solubility} \tab Solubility of calcite (cf. Manning et al., 2013) and \CO2 (cf. Stumm and Morgan, 1996) \cr
    \code{gold} \tab Solubility of gold (Akinfiev and Zotov; 2001; Stef{\aacute}nsson and Seward, 2004; Williams-Jones et al., 2009) \cr
    \code{contour} \tab Gold solubility contours on a log fO2 - pH diagram (Williams-Jones et al., 2009) \cr
    \code{sphalerite} \tab Solubility of sphalerite (Akinfiev and Tagirov, 2014) \cr
    \code{wjd} \tab \eqn{G}{G} minimization: prebiological atmospheres (Dayhoff et al., 1964) and cell periphery of yeast \cr
    \code{dehydration} \tab \logK of dehydration reactions; SVG file contains tooltips and links \cr
    \code{bugstab} \tab Formation potential of microbial proteins in colorectal cancer (Dick, 2016) \cr
    \code{Shh} \tab Affinities of transcription factors relative to Sonic hedgehog (Dick, 2015) \cr
    \code{saturation} \tab Equilibrium activity diagram showing activity ratios and mineral saturation limits (Bowers et al., 1984) \cr
    \code{adenine} \tab HKF regression of heat capacity and volume of aqueous adenine (Lowe et al., 2017) \cr
    \code{DEW} \tab Deep Earth Water (DEW) model for high pressures (Sverjensky et al., 2014a and 2014b) \cr
    \code{lambda} \tab Effects of lambda transition on thermodynamic properties of quartz (Berman, 1988) \cr
    \code{TCA} \tab Standard Gibbs energies of the tricarboxylic (citric) acid cycle (Canovas and Shock, 2016) \cr
    \code{aluminum} \tab Reactions involving Al-bearing minerals (Zimmer et al., 2016; Tutolo et al., 2014) \cr
    \code{carboxylase} \tab Rank abundance distribution for RuBisCO and acetyl-CoA carboxylase \cr
    \code{AkDi} \tab Henry's constant of dissolved gases (Akinfiev and Diamond, 2003) \cr
    \code{comproportionation} \tab Gibbs energy of sulfur comproportionation (Amend et al., 2020) \cr
  }

For either function, if \code{save.png} is TRUE, the plots are saved in \code{\link{png}} files whose names begin with the names of the help topics or demos.

Two of the demos have external dependencies and are not automatically run by \code{demos}.
\samp{dehydration} creates an interactive SVG file; this demo depends on \CRANpkg{RSVGTipsDevice}, which is not available for Windows.
\samp{carboxylase} creates an animated GIF; this demo requires that the ImageMagick \code{convert} commmand be available on the system (tested on Linux and Windows).

\samp{carboxylase} animates diagrams showing rankings of calculated chemical activities along a combined \T and \logaH2 gradient, or makes a single plot on the default device (without conversion to animated GIF) if a single temperature (\code{T}) is specified in the code.
To run this demo, an empty directory named \samp{png} must be present (as a subdirectory of the R working directory).
The proteins in the calculation are 24 carboxylases from a variety of organisms.
There are 12 ribulose phosphate carboxylase and 12 acetyl-coenzyme A carboxylase; 6 of each type are from nominally mesophilic organisms and 6 from nominally thermophilic organisms, shown as blue and red symbols on the diagrams.
The activities of hydrogen at each temperature are calculated using \eqn{\log a_{\mathrm{H_{2}}_{\left(aq\right)}}=-11+3/\left(40\times T\left(^{\circ}C\right)\right)}{logaH2 = -11 + 3/40 * T(degC)}; this equation comes from a model of relative stabilities of proteins in a hot-spring environment (Dick and Shock, 2011).

}

\section{Warning}{
The discontinuities apparent in the plot made by the \code{NaCl} demo illustrate limitations of the "\emph{g} function" for charged species in the revised HKF model (the 355 \degC boundary of region II in Figure 6 of Shock et al., 1992).
Note that \acronym{SUPCRT92} (Johnson et al., 1992) gives similar output at 500 bar.
However, \acronym{SUPCRT} does not output thermodynamic properties above 350 \degC at \Psat; see Warning in \code{\link{subcrt}}.
}

\references{
Akinfiev, N. N. and Diamond, L. W. (2003) Thermodynamic description of aqueous nonelectrolytes at infinite dilution over a wide range of state parameters. \emph{Geochim. Cosmochim. Acta} \bold{67}, 613--629. \url{https://doi.org/10.1016/S0016-7037(02)01141-9}

Akinfiev, N. N. and Tagirov, B. R. (2014) Zn in hydrothermal systems: Thermodynamic description of hydroxide, chloride, and hydrosulfide complexes. \emph{Geochem. Int.} \bold{52}, 197--214. \url{https://doi.org/10.1134/S0016702914030021}

Akinfiev, N. N. and Zotov, A. V. (2001) Thermodynamic description of chloride, hydrosulfide, and hydroxo complexes of Ag(I), Cu(I), and Au(I) at temperatures of 25-500\degC and pressures of 1-2000 bar. \emph{Geochem. Int.} \bold{39}, 990--1006. \url{http://pleiades.online/cgi-perl/search.pl/?type=abstract&name=geochem&number=10&year=1&page=990}

Aksu, S. and Doyle, F. M. (2001) Electrochemistry of copper in aqueous glycine solutions. \emph{J. Electrochem. Soc.} \bold{148}, B51--B57. \url{https://doi.org/10.1149/1.1344532}

Amend, J. P. and Shock, E. L. (1998) Energetics of amino acid synthesis in hydrothermal ecosystems. \emph{Science} \bold{281}, 1659--1662. \url{https://doi.org/10.1126/science.281.5383.1659}

Amend, J. P. and Shock, E. L. (2001) Energetics of overall metabolic reactions of thermophilic and hyperthermophilic Archaea and Bacteria. \emph{FEMS Microbiol. Rev.} \bold{25}, 175--243. \url{https://doi.org/10.1016/S0168-6445(00)00062-0}

Amend, J. P., Aronson, H. S., Macalady, J. and LaRowe, D. E. (2020) Another chemolithotrophic metabolism missing in nature: sulfur comproportionation. \emph{Environ. Microbiol.} \url{https://doi.org/10.1111/1462-2920.14982}

Azadi, M. R., Karrech, A., Attar, M. and Elchalakani, M. (2019) Data analysis and estimation of thermodynamic properties of aqueous monovalent metal-glycinate complexes. \emph{Fluid Phase Equilib.} \bold{480}, 25-40. \url{https://doi.org/10.1016/j.fluid.2018.10.002}

Berman, R. G. (1988) Internally-consistent thermodynamic data for minerals in the system Na{\s2}O-K{\s2}O-CaO-MgO-FeO-Fe{\s2}O{\s3}-Al{\s2}O{\s3}-SiO{\s2}-TiO{\s2}-H{\s2}O-CO{\s2}. \emph{J. Petrol.} \bold{29}, 445-522. \url{https://doi.org/10.1093/petrology/29.2.445}

Bowers, T. S., Jackson, K. J. and Helgeson, H. C. (1984) \emph{Equilibrium Activity Diagrams for Coexisting Minerals and Aqueous Solutions at Pressures and Temperatures to 5 kb and 600\degC}, Springer-Verlag, Berlin, 397 p. \url{http://www.worldcat.org/oclc/11133620}

Canovas, P. A., III and Shock, E. L. (2016) Geobiochemistry of metabolism: Standard state thermodynamic properties of the citric acid cycle. \emph{Geochim. Cosmochim. Acta} \bold{195}, 293--322. \url{https://doi.org/10.1016/j.gca.2016.08.028}

Dayhoff, M. O. and Lippincott, E. R. and Eck, R. V. (1964) Thermodynamic Equilibria In Prebiological Atmospheres. \emph{Science} \bold{146}, 1461--1464. \url{https://doi.org/10.1126/science.146.3650.1461}

Dick, J. M. and Shock, E. L. (2011) Calculation of the relative chemical stabilities of proteins as a function of temperature and redox chemistry in a hot spring. \emph{PLoS ONE} \bold{6}, e22782. \url{https://doi.org/10.1371/journal.pone.0022782}

Dick, J. M. (2015) Chemical integration of proteins in signaling and development. \emph{bioRxiv}. \url{https://doi.org/10.1101/015826}

Dick, J. M. (2016) Proteomic indicators of oxidation and hydration state in colorectal cancer. \emph{PeerJ} \bold{4}:e2238. \url{https://doi.org/10.7717/peerj.2238}

Garrels, R. M. and Christ, C. L. (1965) \emph{Solutions, Minerals, and Equilibria}, Harper & Row, New York, 450 p. \url{http://www.worldcat.org/oclc/517586}

Johnson, J. W., Oelkers, E. H. and Helgeson, H. C. (1992) SUPCRT92: A software package for calculating the standard molal thermodynamic properties of minerals, gases, aqueous species, and reactions from 1 to 5000 bar and 0 to 1000\degC. \emph{Comp. Geosci.} \bold{18}, 899--947. \url{https://doi.org/10.1016/0098-3004(92)90029-Q}

Lowe, A. R., Cox, J. S. and Tremaine, P. R. (2017) Thermodynamics of aqueous adenine: Standard partial molar volumes and heat capacities of adenine, adeninium chloride, and sodium adeninate from \emph{T} = 278.15 K to 393.15 K. \emph{J. Chem. Thermodyn.} \bold{112}, 129--145. \url{https://doi.org/10.1016/j.jct.2017.04.005}

Lu, P. and Zhu, C. (2011) Arsenic Eh--pH diagrams at 25\degC and 1 bar. \emph{Environ. Earth Sci.} \bold{62}, 1673--1683. \url{https://doi.org/10.1007/s12665-010-0652-x}
  
Manning, C. E., Shock, E. L. and Sverjensky, D. A. (2013) The chemistry of carbon in aqueous fluids at crustal and upper-mantle conditions: Experimental and theoretical constraints. \emph{Rev. Mineral. Geochem.} \bold{75}, 109--148. \url{https://doi.org/10.2138/rmg.2013.75.5}

Schulte, M. D. and Shock, E. L. (1995) Thermodynamics of Strecker synthesis in hydrothermal systems. \emph{Orig. Life Evol. Biosph.} \bold{25}, 161--173. \url{https://doi.org/10.1007/BF01581580}

Shock, E. L. and Koretsky, C. M. (1995) Metal-organic complexes in geochemical processes: Estimation of standard partial molal thermodynamic properties of aqueous complexes between metal cations and monovalent organic acid ligands at high pressures and temperatures. \emph{Geochim. Cosmochim. Acta} \bold{59}, 1497--1532. \url{https://doi.org/10.1016/0016-7037(95)00058-8}

Shock, E. L., Oelkers, E. H., Johnson, J. W., Sverjensky, D. A. and Helgeson, H. C. (1992) Calculation of the thermodynamic properties of aqueous species at high pressures and temperatures: Effective electrostatic radii, dissociation constants and standard partial molal properties to 1000 \degC and 5 kbar. \emph{J. Chem. Soc. Faraday Trans.} \bold{88}, 803--826. \url{https://doi.org/10.1039/FT9928800803}

Stef{\aacute}nsson, A. and Seward, T. M. (2004) Gold(I) complexing in aqueous sulphide solutions to 500\degC at 500 bar. \emph{Geochim. Cosmochim. Acta} \bold{68}, 4121--4143. \url{https://doi.org/10.1016/j.gca.2004.04.006}

Stumm, W. and Morgan, J. J. (1996) \emph{Aquatic Chemistry: Chemical Equilibria and Rates in Natural Waters}, John Wiley & Sons, New York, 1040 p. \url{http://www.worldcat.org/oclc/31754493}

Sverjensky, D. A., Harrison, B. and Azzolini, D. (2014a) Water in the deep Earth: The dielectric constant and the solubilities of quartz and corundum to 60 kb and 1,200 \degC. \emph{Geochim. Cosmochim. Acta} \bold{129}, 125--145. \url{https://doi.org/10.1016/j.gca.2013.12.019}

Sverjensky, D. A., Hemley, J. J. and D'Angelo, W. M. (1991) Thermodynamic assessment of hydrothermal alkali feldspar-mica-aluminosilicate equilibria. \emph{Geochim. Cosmochim. Acta} \bold{55}, 989-1004. \url{https://doi.org/10.1016/0016-7037(91)90157-Z}

Sverjensky, D. A., Stagno, V. and Huang, F. (2014b) Important role for organic carbon in subduction-zone fluids in the deep carbon cycle. \emph{Nat. Geosci.} \bold{7}, 909--913. \url{https://doi.org/10.1038/ngeo2291}

Tutolo, B. M., Kong, X.-Z., Seyfried, W. E., Jr. and Saar, M. O. (2014) Internal consistency in aqueous geochemical data revisited: Applications to the aluminum system. \emph{Geochim. Cosmochim. Acta} \bold{133}, 216--234. \url{https://doi.org/10.1016/j.gca.2014.02.036}

Williams-Jones, A. E., Bowell, R. J. and Migdisov, A. A. (2009) Gold in solution. \emph{Elements} \bold{5}, 281--287. \url{https://doi.org/10.2113/gselements.5.5.281}

Zimmer, K., Zhang, Y., Lu, P., Chen, Y., Zhang, G., Dalkilic, M. and Zhu, C. (2016) SUPCRTBL: A revised and extended thermodynamic dataset and software package of SUPCRT92. \emph{Comp. Geosci.} \bold{90}, 97--111. \url{https://doi.org/10.1016/j.cageo.2016.02.013}
}


\examples{
\dontshow{reset()}
\dontshow{opar <- par(no.readonly=TRUE)}
demos(c("ORP", "NaCl"))
\dontshow{par(opar)}
}

\concept{Other tools}
