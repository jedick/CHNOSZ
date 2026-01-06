\name{NEWS}
\title{News for Package 'CHNOSZ'}
\encoding{UTF-8}

% Macros
\newcommand{\degC}{\ifelse{latex}{\eqn{^{\circ}}C}{\ifelse{html}{\out{&deg;}C}{°C}}}
\newcommand{\H2O}{\ifelse{latex}{\eqn{\mathrm{H_{2}O}}}{\ifelse{html}{\out{H<sub>2</sub>O}}{H2O}}}
\newcommand{\Hplus}{\ifelse{latex}{\eqn{\mathrm{H^{+}}}}{\ifelse{html}{\out{H<sup>+</sup>}}{H+}}}
% Subscript and superscript
\newcommand{\s}{\ifelse{latex}{\eqn{_{#1}}}{\ifelse{html}{\out{<sub>#1</sub>}}{#1}}}
\newcommand{\S}{\ifelse{latex}{\eqn{^{#1}}}{\ifelse{html}{\out{<sup>#1</sup>}}{^#1}}}
% Links to vignettes 20220723
\newcommand{\viglink}{\ifelse{html}{\out{<a href="../CHNOSZ/doc/#1.html"><strong>#1.Rmd</strong></a>}}{\bold{#1.Rmd}}}
% Heat capacity and standard Gibbs energy 20231115
\newcommand{\Cp}{\ifelse{latex}{\eqn{C_P}}{\ifelse{html}{\out{<I>C<sub>P</sub></I>}}{Cp}}}
\newcommand{\DG0}{\ifelse{latex}{\eqn{{\Delta}G^{\circ}}}{\ifelse{html}{\out{&Delta;<I>G</I>&deg;}}{ΔG°}}}

\section{Changes in CHNOSZ version 2.2.0-17 (2026-01-06)}{

    \itemize{

      \item \code{label.plot()} and \code{label.figure()} now handle
      logarithmic scales for x and y axes.

      \item Fix bug preventing \code{mosaic()} from being used to form \H2O.
      
      \item For IAPWS95 water model, add 0.1 bar (0.01 MPa) to values of Psat >
      1 bar computed from \code{WP02.auxiliary()} to prevent Born functions,
      which are calculated as finite difference derivatives, from exploding.

      \item Add stibnite (Sb\s{2}S\s{3}) from
      \href{https://doi.org/10.3133/b2131}{Robie and Hemingway (1995)}.

      \item Add \code{phosphorylate()} to calculate affinity of phosphorylation
      reactions with pH-dependent speciation, including Mg-complexed species.

      \item Add \code{hyphen.in.pdf()} for replacing minus signs with
      hyphens in text added to PDF plots.

      \item Use \href{https://online.ucpress.edu/elementa}{Elementa} reference
      style (like \href{https://www.csemanual.org}{CSE} name-year style but all
      authors are listed) in Rd files.

      \item Add colors to \code{demo/NaCl.R}.

      \item Add aqueous SO\s{2} and the trisulfur radical ion (S\s{3}\S{-}) to
      basis species for mosaic calculations in examples, demos, and vignettes.
      However, some lines with S\s{3}\S{-} are commented because of large
      slowdowns (especially in anintro.Rmd and multi-metal.Rmd). Thanks to
      Evgeniy Bastrakov for suggesting this change.

      \item Remove \code{lex()}, which only wrapped its argument in
      \code{as.expression(c())}.

    }

}


\section{Changes in CHNOSZ version 2.2.0 (2025-06-20)}{

  \subsection{OBIGT DATABASE FORMAT}{
    \itemize{

      \item Mineral species whose \code{T} value is a Cp equation limit or
      phase transition temperature are indicated by \samp{CGL} or
      \samp{CGL_Ttr}, respectively, in \code{model}. At temperatures above the
      \code{T} limit, \code{subcrt()} extrapolates values and produces a
      warning for \samp{CGL} species or returns NA values for \samp{CGL_Ttr}
      species. The addition of \samp{CGL_Ttr} replaces the previous strategy of
      using signed values of \code{T}.

    }
  }

  \subsection{OBIGT DEFAULT DATA}{
    \itemize{

      \item Add tantalum species from
      \href{https://doi.org/10.1016/j.gca.2024.10.019}{Hu et al. (2024)}.

      \item Add germanium, indium, gallium, and cadmium minerals and aqueous
      species from the compilation of
      \href{https://doi.org/10.1080/08120099.2024.2423772}{Huston and Bastrakov
        (2024)}. Move Cd and In chloride, In fluoride, and Ga hydroxide
      complexes from \href{https://doi.org/10.1016/S0016-7037(96)00339-0}{Shock
        et al. (1997)} and
      \href{https://doi.org/10.1016/S0016-7037(97)00009-4}{Sverjensky et al.
        (1997)} to optional data. Thanks to Evgeniy Bastrakov for providing
      data files for these species and for Co, Ta, Nb, U, and REE minerals and
      uranyl complexes mentioned below.

      \item Add palladium and platinum hydrosulfide complexes
      (Pd(HS)\s{4}\S{-2} and Pt(HS)\s{4}\S{-2}) from
      \href{https://doi.org/10.1016/j.gca.2024.12.019}{Laskar et al. (2025)}.

      \item Add aqueous uranyl complexes from
      \href{https://doi.org/10.1016/j.gca.2024.04.023}{Migdisov et al. (2025)}.
      Move aqueous uranium species from
      \href{https://doi.org/10.1016/S0016-7037(97)00240-8}{Shock, Sassani and
        Betz (1997)} (except for HUO\s{4}\S{-}) to optional data.

      \item Add alkylamines, benzylamines, and aminiums from
      \href{https://doi.org/10.1016/j.gca.2024.03.013}{Robinson et al. (2024)}.

      \item Add tin chloride complexes and updated parameters for Sn\S{+2} and
      Th\S{+4} from \href{https://doi.org/10.1016/j.gsf.2023.101624}{Liu et al.
        (2023)}.

      \item Add uranium, uraninite (UO\s{2}), beta-UO\s{2.25}, UO\s{2.25},
      alpha-UO\s{2.3333}, beta-UO\s{2.3333}, UO\s{2.6667}, gamma-UO\s{3},
      rutherfordine (UO\s{2}CO\s{3}), beta-UO\s{2}(OH)\s{2}, and
      Na\s{2}U\s{2}O\s{7} from
      \href{https://doi.org/10.1787/bf86a907-en}{Grenthe et al. (2020)}.

      \item Add Fe-bearing serpentinization-relevant minerals (ferrous
      hydroxide, ferrotremolite, greenalite, minnesotaite, Mg-cronstedtite, and
      hisingerite) from \href{https://doi.org/10.2475/03.2020.01}{Leong and
        Shock (2020)}.

      \item Add REE minerals from the compilation of
      \href{https://doi.org/10.1016/j.chemgeo.2016.06.005}{Migdisov et al.
        (2016)} with modifications.
      
      \item Add yttrium fluoride complexes (YF\S{+2} and YF\s{2}\S{+}) from
      \href{https://doi.org/10.1016/j.gca.2013.07.031}{Loges et al. (2013)}.

      \item Add REE fluoride and chloride complexes from
      \href{https://doi.org/10.1016/j.gca.2009.08.023}{Migdisov et al. (2009)}
      (previously available in CHNOSZ 0.9-5 to 1.0.8). Move aqueous REE species
      from \href{https://doi.org/10.1016/0016-7037(95)00314-P}{Haas et al.
        (1995)} to optional data.

      \item Add beta-Ta\s{2}O\s{5} from
      \href{https://doi.org/10.1016/j.jct.2008.12.006}{Jacob et al. (2009)}.

      \item Add Co(OH)\s{2} from
      \href{https://pid.geoscience.gov.au/dataset/ga/25473}{Shvarov and
        Bastrakov, 1999}.

      \item Add chromite (FeCr\s{2}O\s{4}) and cassiterite (SnO\s{2}) from
      \href{https://doi.org/10.3133/b2131}{Robie and Hemingway (1995)}. Move
      cassiterite from
      \href{https://doi.org/10.2113/gsecongeo.80.5.1365}{Jackson and Helgeson
        (1985)} to optional data.

      \item Add spherocobaltite (CoCO\s{3}) from
      \href{https://doi.org/10.1007/978-3-662-02293-1}{Barin et al. (1977)}.

      \item Add CoSO\s{4} and beta-Nb\s{2}O\s{5} from
      \href{https://www.worldcat.org/oclc/8139446}{Naumov et al. (1974)}. A PDF
      of this publication is available from the
      \href{https://ntrl.ntis.gov/NTRL/}{National Technical Reports Library}.
      
      \item Move americium complexes from optional to default data. Entropy of
      the element Am has been listed in \file{element.csv} since CHNOSZ 1.4.0.

      \item Cite \href{https://doi.org/10.1017/CBO9781139167604}{Denbigh
        (1981)} for the hydrogen ion convention.

      \item Due to an apparent typo in
      \href{https://srd.nist.gov/JPCRD/jpcrdS2Vol11.pdf}{Wagman et al. (1982)},
      Δ\emph{H}°\s{f} for In\S{+3} and derived species was corrected by -7335
      cal/mol to maintain GHS consistency. This is in addition to the
      correction for S° previously noted by
      \href{https://doi.org/10.1016/0016-7037(88)90181-0}{Shock and Helgeson
        (1988)}.

      \item Fix sign of \emph{c} parameter for nesquehonite. Thanks to James
      Leong and Grayson Boyer.

      \item Restore steam (H\s{2}O gas).

    }
  }

  \subsection{OBIGT OPTIONAL DATA}{
    \itemize{

      \item Make new directory \file{OBIGT/testing} to hold optional data files
      for testing, including \file{SiO2.csv} (renamed from \file{AS04.csv}) and
      \file{GEMSFIT.csv} in addition to those listed below.

      \item Add \file{testing/IGEM.csv}: Au and Cu complexes and
      SO\s{4}\S{-2}-bearing species from researchers at IGEM RAS.

      \item Add \file{testing/Sandia.csv}: Mg(OH)\s{2}-MgCl\s{2}-H\s{2}O phases from
      \href{https://doi.org/10.1016/j.apgeochem.2024.106032}{Knight et al.
        (2024)}.

    }
  }

  \subsection{DOCUMENTATION}{
    \itemize{

      \item Major revision of \viglink{anintro}.

      \item Add \file{demo/MgATP.R}: speciation of ATP with H\S{+} and
      Mg\S{+2}, based on
      \href{https://search.worldcat.org/title/51242181}{Alberty (2003)} (moved
      from \file{anintro.Rmd}).

      \item Update \file{demo/contour.R} to use conditions taken from
      \href{https://doi.org/10.1016/j.oregeorev.2022.105231}{Ding et al.
        (2023)} and optional data for Au species from IGEM RAS, including
      polynuclear Au\s{2}S\s{2}\S{-2}.

      \item Add \file{demo/rubisco_Zc.R}, moved from \file{anintro.Rmd}.

      \item Add \file{demo/sum_S.R}: log \emph{f}O\s{2}--log \emph{m} ΣS
      diagram for Fe-S-O-C minerals with solubility contours for Fe and Au.

      \item Add \file{demo/uranyl.R}: total (carbonate|sulfate)-pH diagrams for
      uranyl species, after
      \href{https://doi.org/10.1016/j.gca.2024.04.023}{Migdisov et al. (2025)}.

      \item Add FAQ questions:
     
      \itemize{
        
        \item Why are mineral stability boundaries curved on mosaic diagrams?

        \item Where do the names CHNOSZ, OBIGT, and subcrt come from?

        \item What are the main limitations of CHNOSZ?

      }

    }
  }

  \subsection{NEW FEATURES}{
    \itemize{

      \item \code{axis.label()} now uses summation symbol (Σ) for total
      activity or molality of basis species used as a plotting variable in
      \code{mosaic()}.

      \item Restore \code{lines} to the output of \code{diagram()} with the the
      x and y coordinates of stability field boundaries.

      \item \code{check.EOS()} now uses values of Born coefficients \emph{X}
      and \emph{Q} that are consistent with either SUPCRT92 or DEW, depending
      on the \code{model} listed for the species in OBIGT (HKF or DEW,
      respectively).

      \item \code{rank.affinity()} now rescales average ranks of groups to the
      same bounds.

    }
  }

  \subsection{CHANGES TO NAMES AND DEFAULTS}{
    \itemize{

      \item In \code{diagram()}, the default for the \strong{type} argument
      when using the output from \code{solubility()} is now \code{loga.balance}
      (sum of activities of aqueous species) rather than \code{loga.equil}
      (activities of individual aqueous species).

      \item In \code{NaCl()}, rename argument \strong{m_tot} to \strong{m_NaCl}
      (moles of of NaCl added to 1 kg H\s{2}O) and rename output with
      \samp{m_Naplus}, \samp{m_Clminus}, \samp{m_NaCl0}, and \samp{m_HCl0}
      (molalities of aqueous species).

      \item Rename \code{logB.to.OBIGT()} to \code{logK.to.OBIGT()} and add a
      \strong{state} argument which accepts \samp{aq} (aqueous species using
      heat capacity parameters in the \samp{HKF} model; the default) or
      \samp{cr} (minerals using the Maier-Kelley equation in the \samp{CGL}
      model; a new option).

      \item In \code{water()}, rename \strong{P1} logical argument to
      \strong{Psat_floor} numeric argument.

    }
  }

  \subsection{BUG FIXES}{
    \itemize{

      \item \code{water()} with the IAPWS95 model now respects the value of
      \strong{Psat_floor}, making it possible to remove the floor of 1 bar.

      \item Fix using Eh in transect mode of \code{affinity()} (feature was not
      working since at least version 0.9-7).

    }
  }

  \subsection{OTHER CHANGES}{
    \itemize{

      \item \code{expr.species()} no longer subscripts formulas after a
      variable (activity, fugacity, or molality: \emph{a}, \emph{f}, or
      \emph{m}). This avoids readability issues with double subscripting of
      elements' coefficients. An example of the new formatting is \dQuote{log
        \emph{f}O\s{2}}. 
      
      \item Add scripts and data files in \code{extdata/protein/Cas} for amino
      acid compositions of CRISPR-associated (Cas) proteins in different
      classes and subtypes listed by
      \href{https://doi.org/10.1038/s41579-019-0299-x}{Makarova et al. (2020)}.


      \item Remove \code{ispecies} from output of \code{check.OBIGT()} to avoid
      trivial diffs.

    }
  }

  \subsection{REMOVED FEATURES}{
    \itemize{

      \item Remove \code{add.alpha()}. Instead of this function,
      \code{grDevices::adjustcolor()} is now used in \code{stack_mosaic()} to
      add transparency.

      \item Remove \code{seq2aa()}.

      \item Move \code{read.fasta()}, \code{count.aa()}, and \code{aasum()} to
      \href{https://cran.r-project.org/package=canprot}{canprot} with different
      names (\code{read_fasta()}, \code{count_aa()}, and \code{sum_aa()}).

      \item Move \file{equilibrium.Rmd} and \code{protein.equil()} to
      \href{https://github.com/jedick/JMDplots}{JMDplots}.

    }
  }

}

\section{Changes in CHNOSZ version 2.1.0 (2024-02-11)}{

  \subsection{NEW FEATURES}{
    \itemize{

      \item Add \viglink{FAQ} vignette with the following questions: How is
      \sQuote{CHNOSZ} pronounced?, How should CHNOSZ be cited?, What
      thermodynamic models are used in CHNOSZ?, When and why do equal-activity
      boundaries depend on total activity?, How can minerals with polymorphic
      transitions be added to the database?, How can I make a diagram with the
      trisulfur radical ion (S\s{3}\S{-})?, In OBIGT, what is the meaning of
      \code{T} for solids, liquids, and gases?, and How can mineral pH buffers
      be plotted?

      \item Add \strong{use.polymorphs} argument to \code{subcrt()} to allow
      turning off automatic identification of stable polymorphs. This is used
      in the FAQ (How can minerals with polymorphic transitions be added to the
      database?).

      \item Restore \file{EOSregress.R} and its reverse dependencies
      \viglink{eos-regress} and \file{demo/adenine.R}. Thanks to Kris Fecteau
      and Kirt Robinson for the suggestion.

      \item \code{mosaic()} now handles the situation where the alternate basis
      species include one of the axis variables on a diagram (by changing the
      argument names in its internal calls to \code{affinity()}) and adjusts
      the labels for the diagram accordingly (e.g. \dQuote{total C},
      \dQuote{total S}, etc.). Thanks to Evgeniy Bastrakov for the suggestion.

      \item The environment variable CHNOSZ_BUILD_LARGE_VIGNETTES is used to
      control \strong{dpi} in \pkg{knitr} chunk options. Setting this variable
      results in larger vignettes that are used for the CHNOSZ website; if this
      variable is unset (as in CRAN checks), a smaller package is built.

      \item \code{rank.affinity()} now returns average group rankings as
      percentage values.

    }
  }


  \subsection{OBIGT DATABASE}{
    \itemize{

      \item \code{subcrt()} now outputs \DG0 and issues a warning above the
      temperature limit for \Cp equations for minerals and other \samp{CGL}
      species. The temperature limit for \Cp equations is stored as the
      opposite (negative) value in OBIGT. The previous behavior of outputting
      NA for \DG0 and issuing a message is now reserved for phase stability
      limits (i.e., melting, vaporization, or decomposition). See \viglink{FAQ}
      for details.

      \item Add carrollite (Cu\s{0.92}Co\s{2.07}S\s{4}) from
      \href{https://doi.org/10.3749/2300028}{von der Heyden et al. (2024)}.

      \item Modify cattierite (CoS\s{2}) and linnaeite (Co\s{3}S\s{4}) to use
      parameters from \href{https://www.worldcat.org/oclc/16131757}{Pankratz et
        al. (1987)} as listed by
      \href{https://doi.org/10.5382/econgeo.4895}{Williams-Jones and Vasyukova
        (2022)}.
      
      \item Add 36 charged amino acids from
      \href{https://doi.org/10.1007/s11084-015-9428-3}{Kitadai (2015)} and
      update several others (Arg\S{+}, Arg\S{-}, Asp\S{-}, Glu\S{-}, His\S{+},
      Lys\S{+}, Lys\S{-}, and Tyr\S{-}). Thanks to Grayson Boyer.

      \item In \code{OBIGT()}, load \file{Berman_cr.csv} before
      \file{inorganic_cr.csv} so that \code{info("Fe2O3")} matches hematite
      rather than maghemite.

      \item Rename \code{checkGHS()} and \code{checkEOS()} to
      \code{check.GHS()} and \code{check.EOS()} and change the default for
      \strong{return.difference} to TRUE.

    }
  }


  \subsection{BUG FIXES}{
    \itemize{

      \item Fix bugs in \code{subcrt()} where \strong{exceed.Ttr} and
      \strong{IS} arguments were not applied to automatically balanced
      reactions.

      \item Fix bug in \code{dPdTtr()} where OBIGT database parameters were not
      converted to Joules. This allows re-activating the tests for quartz
      (which depend on the computed transition temperature) in
      \file{test-subcrt.R}.

    }
  }

}


\section{Changes in CHNOSZ version 2.0.0 (2023-03-13)}{

  \subsection{MAJOR USER-VISIBLE CHANGES}{
    \itemize{

      \item Units of Joules instead of calories are now used by default for the
      thermodynamic properties output by \code{subcrt()}. That is,
      \code{E.units("J")} is the default setting. Scripts that implicitly
      depend on the previous default setting of \code{E.units("cal")} may need
      to be modified to produce expected output.

      \item A \code{model} column has been added to the OBIGT thermodynamic
      database. This specifies the model for calculating standard thermodynamic
      properties. Currently available models are \samp{H2O} (water), \samp{HKF}
      (revised Helgeson-Kirkham-Flowers model for aqueous species), \samp{CGL}
      (heat capacity equation for crystalline, gas, and liquid species),
      \samp{Berman} (Berman model for minerals), \samp{AD} (Akinfiev-Diamond
      model for aqueous species), and \samp{DEW} (Deep Earth Water model, a
      variation of HKF). All OBIGT data files that can be read using
      \code{add.OBIGT()} must have this column; there is no back-compatibility
      support for data files in the old format.

      \item Backward compatibility for OBIGT data files without an
      \code{E_units} column has been removed.

      \item Order-of-magnitude (OOM) scaling of heat capacity coefficients for
      the CGL model in OBIGT data files has been removed. OOM scaling of HKF
      and DEW model parameters is still in place.

      \item Change license from GPL (>= 2) to GPL-3.

    }
  }

  \subsection{OBIGT DATABASE}{
    \itemize{

      \item Move aqueous CH\s{4} back to \file{organic_aq.csv} and move gaseous
      CH\s{4} from \file{inorganic_gas.csv} to \file{organic_gas.csv}.
      
      \item Modify \code{OBIGT(no.organics = TRUE)} to prune the database of
      organic compounds \emph{except} CH\s{4}. The ability to remove hundreds
      of more complex organic compounds while keeping CH\s{4} is useful for
      some calculations for abiotic reactions in hydrothermal systems.

      \item Merge \file{biotic_aq.csv} into \file{organic_aq.csv}.

      \item \samp{OBIGT/inorganic_cr.csv}: Add aluminum (Al), cobalt monoxide
      (CoO), guite (Co\s{3}O\s{4}), and zinc from
      \href{https://doi.org/10.3133/b2131}{Robie and Hemingway (1995)}.

      \item \samp{OBIGT/inorganic_cr.csv}: Modify linnaeite (Co\s{3}S\s{4}) to
      use 25 \degC parameters from \href{https://doi.org/10.3133/b2131}{Robie
        and Hemingway (1995)} and heat capacity coefficients from
      \href{https://www.worldcat.org/oclc/16131757}{Pankratz et al. (1987)}.

      \item \samp{OBIGT/inorganic_cr.csv}: Modify Co-pentlandite
      (Co\s{9}S\s{8}) to use parameters from
      \href{https://www.worldcat.org/oclc/16131757}{Pankratz et al. (1987)} as
      listed by \href{https://doi.org/10.5382/econgeo.4895}{Williams-Jones and
        Vasyukova (2022)}.

      \item \samp{OBIGT/organic_aq.csv}: Correct formulas of Na(Mal)-, K(Mal)-,
      and Dy(Mal)+. Thanks to Grayson Boyer.

    }
  }

  \subsection{NEW FEATURES: THERMODYNAMIC DATA}{
    \itemize{

      \item Add \code{logB.to.OBIGT()} to fit selected thermodynamic parameters
      (\samp{G} and \samp{S} at 25 \degC and \samp{c1}, \samp{c2}, and
      \samp{omega} HKF coefficients) to formation constants of aqueous species
      as a function of temperature.

      \item Add vignette \viglink{custom_data} with description of database
      format, data-entry conventions, and examples of customizing the
      thermodynamic database using \code{add.OBIGT()}, \code{mod.OBIGT()}, and
      \code{logB.to.OBIGT()}.

    }
  }

  \subsection{NEW FEATURES: DIAGRAMS}{
    \itemize{

      \item Add \code{stack_mosaic()} to create stacked mosaic diagrams, where
      the species formed in the first layer become the basis species for the
      species formed in the second layer. This function implements the
      procedure for making mosaic diagrams for bimetallic systems described in
      \viglink{multi-metal}.

      \item Add \code{loga_aq} argument to \code{mosiac()} (also present in new
      function \code{stack_mosaic()}) to control the activity of mosaiced
      aqueous basis species. This obviates a workaround that was previously
      used in the Mosaic Stacking 2 section of \strong{multi-metal.Rmd}.

      \item Add \code{add.alpha()} to add alpha value to a color specified in
      any notation accepted by \code{col2rgb()}.

    }
  }

  \subsection{NEW FEATURES: PROTEINS}{
    \itemize{

      \item Add function \code{rank.affinity()} to calculate means of affinity
      rankings for specified groups of proteins, intended for evolutionary
      comparisons. The output of the new function can be used by
      \code{diagram()}.

      \item Add an \strong{as.residue} argument to \code{add.protein()} to
      normalize amino acid compositions by protein length.

      \item Add a check to \code{add.protein()} that all new protein IDs
      (protein name + organism name) are unique.

      \item In \code{affinity()}, print a message to indicate whether
      ionization calculations for proteins are being done.

      \item \code{seq2aa()} allows an empty sequence (zero counts of amino
      acids), useful for intializing a data frame of amino acid composition.

    }
  }

  \subsection{NEW FEATURES: OTHER}{
    \itemize{

      \item \code{NaCl()} has been rewritten to include pH dependence (i.e.
      formation of HCl as well as NaCl) by using \code{affinity()} and
      \code{equilibrate()} instead of algebraic equations.

      \item Ion size parameters used in the \samp{Bdot} equation have been
      moved from the code for \code{nonideal()} to
      \file{extdata/thermo/Bdot_acirc.csv}, which is read into
      \code{thermo()$Bdot_acirc} and can be changed by the user.

      \item The two preceding features are used in a new demo \code{yttrium.R}
      to show speciation of Y-Cl complexes as a function of NaCl concentration,
      pH, \emph{T}, and \emph{P}, after
      \href{https://doi.org/10.1016/j.gca.2020.04.015}{Guan et al. (2020)}.

    }
  }


  \subsection{REMOVED FEATURES}{
    \itemize{

      \item \code{cgl()}, \code{hkf()}, and \code{AD()} are no longer exported.
      \code{subcrt()} should be used for all calculations of thermodynamic
      properties.

      \item Remove \code{EOSregress()} and the associated demo and vignette.

      \item Remove parallel calculations in \code{read.fasta()} and \code{count.aa()}.

      \item Remove \code{revisit()}, \code{findit()}, \code{uniprot.aa()},
      \code{nucleic.formula()}, and \code{nucleic.complement()}.

    }
  }

  \subsection{BUG FIXES}{
    \itemize{

      \item Add \code{stringsAsFactors = FALSE} to a \code{cbind()} call in
      \code{mix()}, needed for \R < 4.0.0. Thanks to Xiaojia Jiang for the bug
      report.

      \item In \strong{multi-metal.Rmd}, work around absence of
      \code{hcl.colors()} in \R < 3.6.0.

      \item \code{diagram()}: Improve handling of length > 1 \code{col},
      \code{lty} and \code{lwd} arguments for predominance diagrams. Because
      field boundaries are drawn for the first species, then that species is
      removed before drawing boundaries for the second species, and so on,
      species for which different-colored or -styled lines are desired should
      be placed at the top of the species list. This change makes possible one
      less call to \code{diagram()} in Mosaic Stacking 1 and 2 in
      \strong{multi-metal.Rmd}, and now the entire chalcopyrite field in Mosaic
      Stacking 2 is bounded by a thick orange line, instead of just the border
      with bornite.

      \item Improve handling of non-integer coefficients in
      \code{expr.species()}. For example, the result for FeS1.33 was previously
      incorrectly formatted as FeS\s{.33} but is now correctly formatted as
      FeS\s{1.33}.

    }
  }

  \subsection{OTHER CHANGES}{
    \itemize{

      \item The \code{\var{reaction}} component of the output of
      \code{subcrt()} contains the \code{model} for each species.

      \item \code{subcrt()} produces an error if an aqueous species with the
      \samp{DEW} \code{model} is requested but the DEW water model isn't
      activated.

      \item Add a \strong{zap} argument to \code{mod.OBIGT()} to clear
      parameters of preexisting species (used by \code{logB.to.OBIGT()}).

      \item In \code{affinity()}, make \code{loga.protein = 0} the default
      (changed from -3).

      \item Add tests \file{stack_mosaic.R} and \file{stack_solubility.R}
      (these create PDF files for visual inspection of results).

      \item \code{axis.label()} produces labels with units delimited by
      parentheses instead of a comma.

      \item \code{seq2aa()} now has the \strong{sequence} argument first and a
      default of NA for \strong{protein} (the protein name).

      \item Remove unused \samp{cutoff} value in \code{thermo()$opt}.

      \item Move \strong{ibasis} argument to first position in
      \code{describe.basis()}.

    }
  }

}

\section{Changes in CHNOSZ version 1.4.3 (2022-02-20)}{

  \subsection{UNIT TESTING}{
    \itemize{

      \item Skip tests on CRAN that fail for M1mac and ATLAS (heat capacity and
      the second derivative of \H2O density in \file{tinytest/test-AD.R}).

    }
  }

}

\section{Changes in CHNOSZ version 1.4.2 (2022-02-08)}{

  \subsection{PLANNED API CHANGE}{
    \itemize{

      \item In the next major release of CHNOSZ (to be version 2.0.0), units of
      Joules instead of calories will be used by default for the thermodynamic
      properties output by \code{subcrt()}. That is, \code{E.units("J")} will
      be the default setting. User scripts that implicitly depend on the
      current default setting of \code{E.units("cal")} will need to be modified
      to produce expected results after this change is made.

    }
  }

  \subsection{THERMODYNAMIC CALCULATIONS}{
    \itemize{

      \item The \code{AD()} function, which provides an implementation of the
      \href{https://doi.org/10.1016/S0016-7037(02)01141-9}{Akinfiev and Diamond
        (2003)} equation of state for aqueous nonelectrolytes, now calculates
      entropy, heat capacity, and volume in addition to the previously
      available Gibbs energy.
      
      \item For minerals with thermodynamic parameters in the equations of
      \href{https://doi.org/10.1093/petrology/29.2.445}{Berman (1988)},
      \code{info()} with a numeric argument now lists values of \code{G},
      \code{H}, \code{S}, \code{Cp}, and \code{V} at 25 \degC and 1 bar (these
      were previously shown as \code{NA}). Thanks to Evgeniy Bastrakov for
      suggesting both of these changes.

    }
  }

  \subsection{OBIGT DATABASE}{
    \itemize{

      \item Add H\s{2}WO\s{4}(aq) and Cp coefficients of scheelite (CaWO\s{4})
      from \href{https://doi.org/10.1016/j.chemgeo.2021.120488}{Liu et al.
        (2021)}.

      \item Add aqueous species in the system Ca-Mg-Na-K-Al-Si-O-H-C-Cl from
      \href{https://doi.org/10.1016/j.gca.2016.04.026}{Miron et al. (2016)} and
      \href{https://doi.org/10.2475/07.2017.01}{Miron et al. (2017)} to
      \file{GEMSFIT.csv}.
      
      \item Add pyrobitumen (C\s{54}H\s{42}) from
      \href{https://doi.org/10.1016/j.gca.2008.03.004}{Helgeson et al.  (2009)}
      to \file{organic_cr.csv}.
      
      \item Move Ar, Xe, CH\s{4}, and CO\s{2} from \file{organic_aq.csv} to
      \file{inorganic_aq.csv}.

      \item Fix formula of CaCl\s{2} in \file{DEW.csv}. Thanks to Grayson Boyer.

      \item Remove \file{OldAA.csv} (superseded thermodynamic parameters for
      amino acids). This file is now available in the
      \href{https://github.com/jedick/JMDplots}{JMDplots} package.

    }
  }

  \subsection{OTHER CHANGES}{
    \itemize{

      \item Tests are now run using the \CRANpkg{tinytest} package.

      \item Add an \strong{add} argument to \code{basis()} to allow adding a
      species to an existing set of basis species.

      \item The \code{AkDi()} function has been renamed to \code{AD()}, and all
      variables and data files likewise use the acronym \acronym{AD}. In
      particular, the Akinfiev-Diamond model is activated for an aqueous
      species by using \code{abbrv = "AD"} in \code{thermo()$OBIGT}.

      \item Names of functions, variables, and data files now use capitalized
      \samp{Berman}, not \samp{berman}.

      \item In \file{vignettes/multi-metal.Rmd} (\dQuote{Diagrams with multiple
        metals}), add a link to the associated paper
      (\href{https://doi.org/10.1016/j.acags.2021.100059}{Dick, 2021}).

    }
  }

  \subsection{REMOVED FUNCTIONS}{
    \itemize{

      \item Remove \code{maxdiff()} and \code{expect_maxdiff()}, previously
      used in tests.

      \item Remove \code{read.blast()}, \code{id.blast()},
      \code{write.blast()}, \code{def2gi()}, and example files in
      \file{extdata/bison}.

      \item Remove \code{eqdata()}, a script for extracting data from EQ3/6
      output files.
      
    }
  }

}

\section{Changes in CHNOSZ version 1.4.1 (2021-04-08)}{

  \subsection{NEW FEATURES}{
    \itemize{

      \item Add \samp{demo/Pourbaix.R} (Eh-pH diagrams with isosolubility
      lines, after \href{https://www.worldcat.org/oclc/563921897}{Pourbaix,
      1974}). This demo depends on the following two changes.

      \item Revise \code{solubility()} to enable calculating overall (i.e.
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

      \item \samp{OBIGT/Berman_cr.csv}: Where possible, use standard
      abbreviations for mineral names (including hematite (Hem) and magnetite
      (Mag)), as listed by \href{https://doi.org/10.2138/am.2010.3371}{Whitney
        and Evans (2010)}.

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
      predominant species extracted from the values actually used to make the
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

      \item Rename \samp{demo/berman.R} to \samp{demo/potassium.R}.

      \item Remove unused limSolve package from Suggests in DESCRIPTION.

      \item Add option \code{thermo()$opt$ionize.aa} to control calculation of
      properties of ionized proteins in \code{affinity()}.

      \item Fix bug where protein ionization calculations with
      \code{ionize.aa()} were incorrectly affected by the setting of
      \code{E.units()}.

      \item The Akinfiev-Diamond model is now activated by setting \code{abbrv
        = "AkDi"} for any aqueous species, instead of setting the charge
      (\code{z}) to NA. This simplifies the more common use case of adding a
      custom aqueous species to OBIGT with mostly NA parameters.

      \item Remove \samp{bases2} argument from \code{mosaic()} (old backward
      compatibility feature).
      
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
