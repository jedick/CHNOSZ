\name{NEWS}
\title{News for Package 'CHNOSZ'}
\encoding{UTF-8}

% macros
\newcommand{\degC}{\ifelse{latex}{\eqn{^{\circ}}C}{\ifelse{html}{\out{&deg;}C}{°C}}}
\newcommand{\H2O}{\ifelse{latex}{\eqn{\mathrm{H_{2}O}}}{\ifelse{html}{\out{H<sub>2</sub>O}}{H2O}}}
\newcommand{\Hplus}{\ifelse{latex}{\eqn{\mathrm{H^{+}}}}{\ifelse{html}{\out{H<sup>+</sup>}}{H+}}}
% subscript and superscript
\newcommand{\s}{\ifelse{latex}{\eqn{_{#1}}}{\ifelse{html}{\out{<sub>#1</sub>}}{#1}}}
\newcommand{\S}{\ifelse{latex}{\eqn{^{#1}}}{\ifelse{html}{\out{<sup>#1</sup>}}{^#1}}}

\section{Changes in CHNOSZ version 1.4.0-31 (2021-03-22)}{

  \subsection{NEW FEATURES}{
    \itemize{

      \item Add \samp{demo/Pourbaix.R} (Eh-pH diagrams with isosolubility
      lines, after \href{https://www.worldcat.org/oclc/563921897}{Pourbaix,
      1974}). This demo depends on the following two changes.

      \item Revamp \code{solubility()} to enable calculating overall (i.e.
      minimum) solubility for multiple minerals. Calculations for multiple
      minerals are also now used in \samp{demo/zinc.R}.

      \item Restore \strong{lty.aq} and \strong{lty.cr} arguments to
      \code{diagram()} to control plotting of aq-aq and cr-cr field boundaries.

      \item Add \strong{no.organics} argument to \code{OBIGT()} to load the
      thermodynamic database without organic species.

    }
  }

  \subsection{DATA UPDATES}{
    \itemize{

      \item \samp{OBIGT/Berman_cr.csv}: Use standard abbreviations for hematite
      (Hem) and magnetite (Mag).

      \item \samp{OBIGT/organic_aq.csv}: Add properties of aqueous fatty acids,
      saccharides, and other species, and add associated \samp{demo/E_coli.R}
      (Gibbs energy of biomass synthesis in \emph{E. coli}, after
      \href{https://doi.org/10.1038/ismej.2015.227}{LaRowe and Amend, 2016}).

      \item \samp{OBIGT/inorganic_aq.csv}: Add properties of aqueous Nb and Ta
      species from \href{https://doi.org/10.1016/j.gca.2020.04.009}{Akinfiev et
        al., 2020}.
      
      \item \samp{OBIGT/inorganic_cr.csv}: Add willemite (Zn\s{2}SiO\s{4}) from
      \href{https://doi.org/10.3133/b2131}{Robie and Hemingway, 1995} and
      \href{https://doi.org/10.1007/978-3-662-02293-1}{Barin et al., 1977}
      (heat capacity equation).

      \item \samp{OBIGT/inorganic_cr.csv}: Move alunite here from
      \samp{SUPCRT92.csv}.

      \item \samp{thermo/protein.csv}: Rename \samp{CDC2_HUMAN} to
      \samp{CDK1_HUMAN} (UniProt:
      \href{https://www.uniprot.org/uniprot/P06493}{P06493}).

    }
  }


  \subsection{OTHER CHANGES}{
    \itemize{

      \item Automatic detection of dissociation reactions was fragile and has
      been removed from \code{solubility()}. The new default (\code{dissociate
        = FALSE}) is to not consider dissociation reactions.

      \item In the return value of \code{diagram()}, \code{predominant.values}
      previously only contained affinities (extracted from the \code{values}
      element of the \code{eout} argument). Now it has the values for the
      predominant species extracted from the actual values used to make the
      plot (the \code{plotvals} list in the output), i.e.  affinities divided
      by the balancing coefficients if \code{eout} is the output of
      \code{affinity}, or activities if \code{eout} is the output of
      \code{equilibrate}. Therefore it can now be used to draw contours or a
      color image showing the activities of the predominant species. This is
      used for a diagram in a new preprint
      (\href{https://doi.org/10.1101/2021.01.29.428804}{Dick, 2021}).

      \item Revise multi-metal.Rmd: Improve mineral abbreviations and placement
      of labels; use updated DFT energies from Materials Project; add Mosaic
      Stacking 2 (minerals and aqueous species); add \emph{K}\S{eff}
      calculation (\href{https://doi.org/10.1016/j.gca.2021.01.038}{Robinson et
      al., 2021}); add Δ\emph{G}\s{pbx} color scale.

      \item Rename \samp{demo/berman.R} to \samp{potassium.R}.

      \item Remove unused limSolve package from Suggests in DESCRIPTION.

      \item Add option \code{thermo()$opt$ionize.aa} to control calculation of
      properties of ionized proteins in \code{affinity()}.

    }
  }

}

\section{Changes in CHNOSZ version 1.4.0 (2020-11-11)}{

  \subsection{MAJOR CHANGES}{
    \itemize{

      \item An \strong{add} argument has been added to \code{species()}.
      \emph{With the default of \code{add = FALSE}, loading new species now
        causes any existing species definition to be \strong{deleted}.} To add
      species to an existing system, use \samp{add = TRUE}. With this change,
      scripts can be made more concise because \code{species(delete =
        TRUE)} is no longer needed to clear the species definition in a series
      of calculations for different systems. 

      \item The naming convention for methane has changed. The name
      \strong{methane} now applies exclusively to the gas; the formula
      \strong{CH4} without a state continues to represent the aqueous species.
      This behavior is consistent with inorganic gases but differs from most
      organic substances, where the name refers to the aqueous species. This
      change was made because in the past there was no way to use \code{info()}
      or \code{subcrt()} with a single character argument to identify gaseous
      methane, which is common in geochemistry. (Note that \code{info("CH4",
        "gas")} can still be used.) A hack to revert to the old behavior, where
      \strong{methane} refers primarily to the aqueous species, is
      \code{mod.OBIGT(info("CH4"), name = "methane")}.

      \item The all-uppercase acronym \strong{OBIGT} is used everywhere it
      appears in file, function, and object names. In particular, the
      thermodynamic database now is located at \code{thermo()$OBIGT}, and
      functions \code{add.OBIGT()} and \code{mod.OBIGT()} replace the previous
      \code{add.obigt()} and \code{mod.obigt()}.

      \item The new default of \code{NULL} for \strong{limit.water} in
      \code{diagram()} causes the main diagram to be overlaid on the fields of
      water stability (white) and instability (gray). Change this to
      \code{TRUE} to obtain the previous default of masking all parts of the
      diagram outside of the water stability field.

    }
  }

  \subsection{NEW FEATURES}{
    \itemize{

      \item Add function \strong{mash()} for combining two diagrams for
      different systems (i.e., simple overlay of diagrams for two single-metal
      systems).

      \item Add function \strong{mix()} for combining two single-metal diagrams
      with a third diagram for bimetallic species. This can be used to produce
      diagrams for a binary system with fixed composition of the metals.

      \item Add function \strong{rebalance()} for making a new diagram by
      secondary balancing between two systems.

      \item Add a \strong{stable} argument to \code{mosaic()} to use previously
      calculated stabilities of species (e.g. minerals) for the changing basis
      species. This allows \code{mosaic()} calculations to be linked in series
      (or "stacked"), for instance to sequentially add metals (Fe, then Cu) to
      a diagram.

      \item Add vignette \strong{multi-metal.Rmd} for examples that use these
      new features to make diagrams for systems with multiple metals
      (specifically Fe-V and Fe-Cu).

    }
  }

  \subsection{CHANGES TO OBIGT DATABASE}{
    \itemize{

      \item Dates in the \samp{date} column have been converted to use ISO 8601
      extended format (e.g. 2020-08-16). The accessory function \code{today()},
      to render the current date using the previous format, has been removed.

      \item Update species parameters for the 2019 version of the Deep Earth
      Water (DEW) model (\href{https://doi.org/10.1016/j.gca.2019.03.027}{Huang
      and Sverjensky, 2019}).

      \item \samp{inorganic_aq.csv}, \samp{inorganic_aq.csv}: Add arsenopyrite,
      loellingite, westerveldite, and aqueous As(OH)\s{3} and AsO(OH)\s{3} from
      \href{https://doi.org/10.1016/j.gca.2007.11.017}{Perfetti et al., 2008}.

      \item \samp{inorganic_aq.csv}: Add Co-chloride species from
      \href{https://doi.org/10.1016/j.gca.2010.12.002}{Liu et al., 2018}.

      \item \samp{inorganic_aq.csv}: Add metal carbonate and bicarbonate
      complexes and FeSO\s{4} from
      \href{https://doi.org/10.1021/acsearthspacechem.9b00016}{St Clair et al.,
        2019}.

      \item \samp{inorganic_cr.csv}: Add wustite, manganese, and Mn oxides
      (manganosite, pyrolusite, bixbyite, and hausmannite) from
      \href{https://doi.org/10.3133/b2131}{Robie and Hemingway,
        1995}.

      \item \samp{inorganic_cr.csv}: Add amorphous Mn(OH)\s{2} from
      \href{https://srd.nist.gov/JPCRD/jpcrdS2Vol11.pdf}{Wagman et al., 1982}
      (GHS) and \href{https://doi.org/10.1016/S0925-8388(98)00739-7}{Senoh et
        al., 1998} (estimated Cp values).

      \item \samp{inorganic_cr.csv}: Add huebnerite (MnWO\s{4}) from
      \href{https://doi.org/10.1016/j.oregeorev.2019.103289}{Liu and Xiao,
        2020}.

      \item \samp{inorganic_cr.csv}: Add linnaeite (Co\s{3}S\s{4}) and
      Co-pentlandite (Co\s{9}S\s{8}) from
      \href{https://www.worldcat.org/oclc/801856147}{Vaughan and Craig, 1978}
      and cattierite (CoS\s{2}) and cobalt from
      \href{https://doi.org/10.3133/b2131}{Robie and Hemingway,
        1995}.

      \item Add GHS and Cp at 25 \degC for almandine, dickite, fluorphlogopite,
      glaucophane, grunerite, halloysite, larnite, pyrope (\samp{SUPCRT92.csv})
      and bromellite (\samp{inorganic_cr.csv}) from
      \href{https://doi.org/10.3133/b2131}{Robie and Hemingway,
        1995}.

      \item \samp{inorganic_cr.csv}: Correct formulas for number of \H2O in
      some As-bearing minerals.

      \item \samp{SUPCRT92.csv}: Change lowercase to uppercase "A"
      (Angstrom) in names of minerals with suffixes indicating interlayer
      spacing (amesite, clinochlore, cronstedtite, daphnite).

      \item Remove bracketed suffixes in references (e.g. [S92] was used to
      indicate that the data first appeared in \acronym{SUPCRT92}).

      \item Use uncompressed data files (i.e. \samp{.csv}, not \samp{.csv.xz})
      and rename \samp{DEW_aq.csv} to \samp{DEW.csv}.

      \item The package is now installed with the bibtex file for the OBIGT
      database (\samp{docs/OBIGT.bib}), which is used by the \pkg{logKcalc}
      package (\url{https://github.com/jedick/logKcalc}).

      \item \samp{thermo/element.csv}: Add entropies for Am, Pu, Cm, Np from
      Thoenen et al., 2014 (\href{https://www.psi.ch/en/les/database}{The
        PSI/Nagra Chemical Thermodynamic Database 12/07}).

    }
  }

  \subsection{DEMOS AND VIGNETTES}{
    \itemize{

      \item Add \samp{demo/comproportionation.R}: Gibbs energy of sulfur
      comproportionation, after
      \href{https://doi.org/10.1111/1462-2920.14982}{Amend et al., 2020}.

      \item Add \samp{demo/zinc.R}: Solubilities of multiple minerals (zincite
      and sphalerite).

      \item Revise and simplify \samp{anintro.Rmd} and \samp{equilibrium.Rmd}
      (converted from \samp{equilibrium.Rnw}).

      \item Revise \samp{OBIGT.Rmd} to reduce the size of the HTML file and
      make deep linking to individual sections work.

      \item Add \samp{viglink} Rd macro so HTML versions of Rd files can link
      to vignettes.

      \item Remove vignette \samp{hotspring.Rnw} and demos \samp{bison.R} and
      \samp{bugstab.R}. The material is now in the \samp{bison.Rmd} vignette
      and \samp{bugstab.R} demo in the \pkg{JMDplots} package
      (\url{https://github.com/jedick/JMDplots}).

      \item New demo \samp{berman.R}, extracted from \samp{berman.Rd}.

      \item Revise demo \samp{mosaic.R} to use \code{mash()} to show S and C
      species together, make better use of transparency, and \code{dy} argument
      to adjust positions of labels.

    }
  }

  \subsection{DIAGRAM IMPROVEMENTS}{
    \itemize{

      \item Change default resolution in \code{affinity()} from 128 to 256.
      This greatly improves the appearance of diagrams made with the default
      settings.

      \item \code{which.pmax()} was rewritten to speed up identification of
      predominant species, and the drawing of field boundaries in \code{diagram()}
      was made more efficient. The latter change has prompted the removal of
      the \samp{lty.aq} and \samp{lty.cr} arguments from \code{diagram()}.

      \item Add a \strong{bottom} argument to \code{ratlab()} to allow changing
      the ion in the denominator to something other than \Hplus.

      \item The \samp{srt}, \strong{dx} (new), and \samp{dy} arguments in
      \code{diagram()} can now be used to rotate and adjust the position of
      field labels, not only line labels. This and other arguments
      (\samp{cex}, \samp{col}, \samp{col.names}, \samp{font},
      \samp{family}, \samp{bold}, \samp{italic}) can have length > 1 to
      apply different settings to each species.

      \item Add a \strong{min.area} argument to \code{diagram()} to specify the
      minimum area of fields that should be labeled. This is useful for
      removing labels from small fields on crowded diagrams.
      
      \item The list returned by \code{diagram()} now includes a
      \strong{predominant.values} component, which has the affinities of the
      predominant species at each grid point. This can be used as shown in the
      \code{multi-metal.Rmd} vignette to compute the Pourbaix energy
      (Δ\emph{G}\s{pbx}) for a metastable material.

    }
  }

  \subsection{OTHER CHANGES}{
    \itemize{

      \item \code{subcrt()}: replace \samp{action.unbalanced} argument with
      \strong{autobalance}, which now provides the ability to prevent
      autobalancing.

      \item Setting the water model with \code{water()} updates the literature
      references in \code{thermo()$OBIGT}.

      \item \code{thermo.refs()} shows CHNOSZ version and date.

      \item \code{subcrt()} and \code{affinity()} use degree symbol
      (\strong{°C}) in messages.

      \item Change \code{thermo$...} to \code{thermo()$...} in messages and
      comments.

      \item \code{mosaic()} now allows a \strong{blend} argument of length > 1 to
      apply a specific setting to each group of basis species.

      \item Fix \samp{palply.Rd} for new warning about \dQuote{Non-file
        package-anchored link(s) in documentation object} in \command{R CMD
        check}.

      \item Remove \code{wjd()} (implementation of Gibbs energy minimization
      method from \href{https://doi.org/10.1063/1.1744264}{White et al., 1958})
      and demo \code{wjd.R} and supporting data file.

      \item Add "QCa" as a keyword for preset species in \code{basis()}
      (glutamine, cysteine, acetic acid, \H2O, O\s{2}).

      \item Add \strong{T} and \strong{P} arguments to \code{retrieve()} to
      require that species have non-NA values of standard Gibbs energy at this
      temperature and pressure.

      \item Change defaults in \code{mod.buffer()} to be more suitable for
      mineral buffers (\code{state = "cr"} and \code{logact = 0}).

      \item Remove extdata/supcrt (scripts to read SUPCRT data files).

      \item Convert this NEWS file to Rd format.

    }
  }

}
