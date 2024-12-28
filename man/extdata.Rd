\encoding{UTF-8}
\name{extdata}
\alias{extdata}
\docType{data}
\title{Extra data}

\description{
  The files in the subdirectories of \code{extdata} provide additional thermodynamic data and other data to support the examples in the package documentation and vignettes. 
  See \code{\link{thermo}} for a description of the files in \code{extdata/OBIGT}, which are used to generate the thermodynamic database.
}

\details{

  Files in \code{Berman} contain thermodynamic data for minerals using the Berman formulation:
  \itemize{
    \item \code{Ber88_1988.csv} contains thermodynamic data for minerals taken from Berman (1988).
    \item Other files with names like \code{xxx_yyyy.csv} contain thermodynamic data from other sources; xxx in the filename corresponds to the reference in \code{\link{thermo}$OBIGT} and yyyy gives the year of publication.
    \code{\link{Berman}} uses these data for the calculation of thermodynamic properties at specified \P and \T, which are then available for use in \code{\link{subcrt}}.
    If there are any duplicated mineral names in the files, only the most recent data are used, as determined by the year in the file name.
    Following conventions used SUPCRT92 (see Helgeson et al., 1978), the names of sanidine and microcline were changed to K-feldspar,high and K-feldspar,low (by using the same names in all data files, loading the optional SUPCRT92 data file updates these minerals rather than makes new ones).
    \item \code{sympy.R} is an R script that uses \CRANpkg{rSymPy} to symbolically integrate Bermans's equations for heat capacity and volume to write experessions for enthalpy, entropy and Gibbs energy.
    \item The \code{testing} directory contains data files based on Berman and Aranovich (1996). These are used to demonstrate the addition of data from a user-supplied file (see \code{\link{Berman}}).
  }

  Files in \code{protein} contain protein sequences and amino acid compositions for proteins.
  \itemize{
    \item \code{rubisco.fasta} Sequences of Rubisco obtained from UniProt (see Dick, 2014). See the vignette \viglink{anintro} for an example that uses this file.
    \item \code{POLG.csv}
      Amino acid compositions of a few proteins used for some tests and examples.
      These are various subunits of the Poliovirus type 1 polyprotein (POLG_POL1M in UniProt).
    \item \code{TBD+05.csv} lists genes with transcriptomic expression changes in carbon limitation stress response experiments in yeast (Tai et al., 2005).
    \item \code{TBD+05_aa.csv} has the amino acid compositions of proteins coded by those genes.
      The last two files are used in \code{demo{"rank.affinity"}}.
  }

  Files in \code{taxonomy} contain taxonomic data files:
  \itemize{
    \item \code{names.dmp} and \code{nodes.dmp} are excerpts of NCBI taxonomy files (\url{https://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz}, accessed 2010-02-15). These files contain only the entries for \emph{Escherichia coli} K-12, \emph{Saccharomyces cerevisiae}, \emph{Homo sapiens}, \emph{Pyrococcus furisosus} and \emph{Methanocaldococcus jannaschii} (taxids 83333, 4932, 9606, 186497, 243232) and the higher-ranking nodes (genus, family, etc.) in the respective lineages. See \code{\link{taxonomy}} for examples that use these files.
  }

  Files in \code{misc} contain miscellaneous data:
  \itemize{
    \item \code{PM90.csv} Heat capacities of four unfolded aqueous proteins taken from Privalov and Makhatadze, 1990. Temperature in \degC is in the first column, and heat capacities of the proteins in J mol\eqn{^{-1}}{^-1} K\eqn{^{-1}}{^-1} in the remaining columns. See \code{\link{ionize.aa}} and the vignette \viglink{anintro} for examples that use this file.
    \item \code{RH95.csv} Heat capacity data for iron taken from Robie and Hemingway, 1995. Temperature in Kelvin is in the first column, heat capacity in J K\eqn{^{-1}}{^-1} mol\eqn{^{-1}}{^-1} in the second. See \code{\link{subcrt}} for an example that uses this file.
    \item \code{SOJSH.csv} Experimental equilibrium constants for the reaction NaCl(aq) = Na+ + Cl- as a function of temperature and pressure taken from Fig. 1 of Shock et al., 1992. See \code{demo("NaCl")} for an example that uses this file.
    \item \code{HWM96_V.csv}, \code{HW97_Cp.csv} Apparent molar volumes and heat capacities of \CH4, \CO2, \H2S, and \NH3 in dilute aqueous solutions reported by Hnědkovský et al., 1996 and Hnědkovský and Wood, 1997. Units are Kelvin, MPa, J/K/mol, and cm3/mol. See \code{demo("AD")}, \code{\link{EOSregress}} and the vignette \viglink{eos-regress} for examples that use these files.
    \item \code{SC10_Rainbow.csv} Values of temperature (\degC), pH and logarithms of activity of \CO2, \H2, \NH4plus, \H2S and \CH4 for mixing of seawater and hydrothermal fluid at Rainbow field (Mid-Atlantic Ridge), taken from Shock and Canovas, 2010. See the vignette \viglink{anintro} for an example that uses this file.
    \item \code{SS98_Fig5a.csv}, \code{SS98_Fig5b.csv} Values of logarithm of fugacity of \O2 and pH as a function of temperature for mixing of seawater and hydrothermal fluid, digitized from Figs. 5a and b of Shock and Schulte, 1998. See the vignette \viglink{anintro} for an example that uses this file.
    \item \code{rubisco.csv} UniProt IDs for Rubisco, ranges of optimal growth temperature of organisms, domain and name of organisms, and URL of reference for growth temperature, from Dick, 2014. See \code{\link{rank.affinity}} and the vignette \viglink{anintro} for examples that use this file.
    \item \code{bluered.txt} Blue - light grey - red color palette, computed using \CRANpkg{colorspace}\code{::diverge_hcl(1000,} \code{c = 100, l = c(50, 90), power = 1)}. This is used by \code{\link{ZC.col}}.
    \item \code{AD03_Fig1?.csv} Experimental data points digitized from Figure 1 of Akinfiev and Diamond, 2003, used in \code{demo("AD")}.
    \item \code{TKSS14_Fig2.csv} Experimental data points digitized from Figure 2 of Tutolo et al., 2014, used in \code{demo("aluminum")}.
    \item \code{Mer75_Table4.csv} Values of log(aK+/aH+) and log(aNa+/aH+) from Table 4 of Merino, 1975, used in \code{demo("aluminum")}.
    \item \code{BZA10.csv} contains supplementary thermodynamic data taken from Bazarkina et al. (2010). The data can be added to the database in the current session using \code{\link{add.OBIGT}}. See \code{\link{add.OBIGT}} for an example that uses this file.
    \item \code{OBIGT_check.csv} contains the results of running \code{\link{check.OBIGT}} to check the internal consistency of entries in the default and optional datafiles.
    \item \code{RH98_Table15.csv} Group stoichiometries for high molecular weight crystalline and liquid organic compounds taken from Table 15 of Richard and Helgeson, 1998. The first three columns have the \code{compound} name, \code{formula} and physical \code{state} (\samp{cr} or \samp{liq}). The remaining columns have the numbers of each group in the compound; the names of the groups (columns) correspond to species in \code{\link{thermo}$OBIGT}. The compound named \samp{5a(H),14a(H)-cholestane} in the paper has been changed to \samp{5a(H),14b(H)-cholestane} here to match the group stoichiometry given in the table. See \code{\link{RH2OBIGT}} for a function that uses this file.
    \item \code{SK95.csv} contains thermodynamic data for alanate, glycinate, and their complexes with metals, taken from Amend and Helgeson (1997) and Shock and Koretsky (1995) as corrected in slop98.dat. These data are used in \code{demo("copper")} and \code{demo("glycinate")}.
    \item \code{LA19_test.csv} contains thermodynamic data for dimethylamine and trimethylamine from LaRowe and Amend (2019) in energy units of both J and cal. This file is used in \code{test-util.data.R}) to check the messages produced by \code{\link{check.GHS}} and \code{\link{check.EOS}}.
  }

}

\references{
Akinfiev, N. N. and Diamond, L. W. (2003) Thermodynamic description of aqueous nonelectrolytes at infinite dilution over a wide range of state parameters. \emph{Geochim. Cosmochim. Acta} \bold{67}, 613--629. \doi{10.1016/S0016-7037(02)01141-9}

Amend, J. P. and Helgeson, H. C. (1997) Calculation of the standard molal thermodynamic properties of aqueous biomolecules at elevated temperatures and pressures. Part 1. L-\alpha-amino acids. \emph{J. Chem. Soc., Faraday Trans.} \bold{93}, 1927--1941. \doi{10.1039/A608126F}

Bazarkina, E. F., Zotov, A. V. and Akinfiev, N. N. (2010) Pressure-dependent stability of cadmium chloride complexes: Potentiometric measurements at 1–1000 bar and 25°C. \emph{Geol. Ore Deposits} \bold{52}, 167--178. \doi{10.1134/S1075701510020054}

Berman, R. G. (1988) Internally-consistent thermodynamic data for minerals in the system Na{\s2}O-K{\s2}O-CaO-MgO-FeO-Fe{\s2}O{\s3}-Al{\s2}O{\s3}-SiO{\s2}-TiO{\s2}-H{\s2}O-CO{\s2}. \emph{J. Petrol.} \bold{29}, 445-522. \doi{10.1093/petrology/29.2.445}

Berman, R. G. and Aranovich, L. Ya. (1996) Optimized standard state and solution properties of minerals. I. Model calibration for olivine, orthopyroxene, cordierite, garnet, and ilmenite in the system FeO-MgO-CaO-Al{\s2}O{\s3}-TiO{\s2}-SiO{\s2}. \emph{Contrib. Mineral. Petrol.} \bold{126}, 1-24. \doi{10.1007/s004100050233}

Dick, J. M. (2014) Average oxidation state of carbon in proteins. \emph{J. R. Soc. Interface} \bold{11}, 20131095. \doi{10.1098/rsif.2013.1095}

Gattiker, A., Michoud, K., Rivoire, C., Auchincloss, A. H., Coudert, E., Lima, T., Kersey, P., Pagni, M., Sigrist, C. J. A., Lachaize, C., Veuthey, A.-L., Gasteiger, E. and Bairoch, A. (2003) Automatic annotation of microbial proteomes in Swiss-Prot. \emph{Comput. Biol. Chem.} \bold{27}, 49--58. \doi{10.1016/S1476-9271(02)00094-4}

Helgeson, H. C., Delany, J. M., Nesbitt, H. W. and Bird, D. K. (1978) Summary and critique of the thermodynamic properties of rock-forming minerals. \emph{Am. J. Sci.} \bold{278-A}, 1--229. \url{https://www.worldcat.org/oclc/13594862}

Hnědkovský, L., Wood, R. H. and Majer, V. (1996) Volumes of aqueous solutions of \CH4, \CO2, \H2S, and \NH3 at temperatures from 298.15 K to 705 K and pressures to 35 MPa. \emph{J. Chem. Thermodyn.} \bold{28}, 125--142. \doi{10.1006/jcht.1996.0011}

Hnědkovský, L. and Wood, R. H. (1997) Apparent molar heat capacities of aqueous solutions of \CH4, \CO2, \H2S, and \NH3 at temperatures from 304 K to 704 K at a pressure of 28 MPa. \emph{J. Chem. Thermodyn.} \bold{29}, 731--747. \doi{10.1006/jcht.1997.0192}

Joint Genome Institute (2007) Bison Pool Environmental Genome. Protein sequence files downloaded from IMG/M (\url{https://img.jgi.doe.gov/})

LaRowe, D. E. and Amend, J. P. (2019) The energetics of fermentation in natural settings. \emph{Geomicrobiol. J.} \bold{36}, 492--505. \doi{10.1080/01490451.2019.1573278}

Merino, E. (1975) Diagenesis in teriary sandstones from Kettleman North Dome, California. II. Interstitial solutions: distribution of aqueous species at 100&deg;C and chemical relation to diagenetic mineralogy. \emph{Geochim. Cosmochim. Acta} \bold{39}, 1629--1645. \doi{10.1016/0016-7037(75)90085-X}

Privalov, P. L. and Makhatadze, G. I. (1990) Heat capacity of proteins. II. Partial molar heat capacity of the unfolded polypeptide chain of proteins: Protein unfolding effects. \emph{J. Mol. Biol.} \bold{213}, 385--391. \doi{10.1016/S0022-2836(05)80198-6}

Richard, L. and Helgeson, H. C. (1998) Calculation of the thermodynamic properties at elevated temperatures and pressures of saturated and aromatic high molecular weight solid and liquid hydrocarbons in kerogen, bitumen, petroleum, and other organic matter of biogeochemical interest. \emph{Geochim. Cosmochim. Acta} \bold{62}, 3591--3636. \doi{10.1016/S0016-7037(97)00345-1}

Robie, R. A. and Hemingway, B. S. (1995) \emph{Thermodynamic Properties of Minerals and Related Substances at 298.15 K and 1 Bar (\eqn{10^5} Pascals) Pressure and at Higher Temperatures}. U. S. Geol. Surv., Bull. 2131, 461 p. \url{https://www.worldcat.org/oclc/32590140}

Shock, E. and Canovas, P. (2010) The potential for abiotic organic synthesis and biosynthesis at seafloor hydrothermal systems. \emph{Geofluids} \bold{10}, 161--192. \doi{10.1111/j.1468-8123.2010.00277.x}

Shock, E. L. and Koretsky, C. M. (1995) Metal-organic complexes in geochemical processes: Estimation of standard partial molal thermodynamic properties of aqueous complexes between metal cations and monovalent organic acid ligands at high pressures and temperatures. \emph{Geochim. Cosmochim. Acta} \bold{59}, 1497--1532. \doi{10.1016/0016-7037(95)00058-8}

Shock, E. L., Oelkers, E. H., Johnson, J. W., Sverjensky, D. A. and Helgeson, H. C. (1992) Calculation of the thermodynamic properties of aqueous species at high pressures and temperatures: Effective electrostatic radii, dissociation constants and standard partial molal properties to 1000 \degC and 5 kbar. \emph{J. Chem. Soc. Faraday Trans.} \bold{88}, 803--826. \doi{10.1039/FT9928800803}

Shock, E. L. and Schulte, M. D. (1998) Organic synthesis during fluid mixing in hydrothermal systems. \emph{J. Geophys. Res.} \bold{103}, 28513--28527. \doi{10.1029/98JE02142}

Tai, S. L., Boer, V. M., Daran-Lapujade, P., Walsh, M. C., de Winde, J. H., Daran, J.-M. and Pronk, J. T. (2005) Two-dimensional transcriptome analysis in chemostat cultures: Combinatorial effects of oxygen availability and macronutrient limitation in \emph{Saccharomyces cerevisiae}. \emph{J. Biol. Chem.} \bold{280}, 437--447. \doi{10.1074/jbc.M410573200}

Tutolo, B. M., Kong, X.-Z., Seyfried, W. E., Jr. and Saar, M. O. (2014) Internal consistency in aqueous geochemical data revisited: Applications to the aluminum system. \emph{Geochim. Cosmochim. Acta} \bold{133}, 216--234. \doi{10.1016/j.gca.2014.02.036}

}

\concept{Thermodynamic data}
