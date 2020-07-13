\name{NEWS}
\title{News for Package 'CHNOSZ'}
\encoding{UTF-8}

% macros
\newcommand{\H2O}{\ifelse{latex}{\eqn{\mathrm{H_{2}O}}}{\ifelse{html}{\out{H<sub>2</sub>O}}{H2O}}}

\section{Changes in CHNOSZ version 1.3.6-31 (2020-07-13)}{

  \subsection{MAJOR CHANGES}{
    \itemize{

      \item The naming convention for methane has changed. The name
      \strong{methane} now applies exclusively to the gas; the formula
      \strong{CH4} without a state continues to represent the aqueous species.
      This behavior is consistent with inorganic gases but differs from most
      organic substances, where the English name refers to the aqueous species.
      This change was made because in the past there was no way to use
      \code{info()} or \code{subcrt()} with a single character argument to
      identify gaseous methane, which is common in geochemistry.  (Note that
      \code{info("CH4", "gas")} can still be used.) A hack to revert to the old
      behavior, where \strong{methane} refers primarily to the aqueous species,
      is \code{mod.OBIGT(info("CH4"), name = "methane")}.

      \item The all-uppercase acronym \strong{OBIGT} is used everywhere it
      appears in file, function, and object names. In particular, the
      thermodynamic database now sits at \code{thermo()$OBIGT}, and functions
      \code{add.OBIGT()} and \code{mod.OBIGT()} replace the previous
      \code{add.obigt()} and \code{mod.obigt()}.

    }
  }

  \subsection{CHANGES TO OBIGT DATABASE}{
    \itemize{

      \item Update species parameters for the 2019 Deep Earth Water (DEW)
      model.

      \item \samp{OBIGT/SUPCRT92.csv}: Change lowercase to uppercase "A"
      (Angstrom) in names of minerals with suffixes indicating interlayer
      spacing (amesite, clinochlore, cronstedtite, daphnite).

      \item \samp{OBIGT/inorganic_aq.csv}: Update arsenopyrite with data from
      Naumov et al., 1974 (via Unitherm).

      \item \samp{OBIGT/inorganic_cr.csv}: Correct syntax for number of \H2O in
      formulas of some As-bearing minerals.

      \item Install bibtex file (\samp{docs/OBIGT.bib}) for use in the
      \pkg{logKcalc} package (\url{https://github.com/jedick/logKcalc}).

      \item Use uncompressed data files (i.e. \samp{.csv}, not \samp{.csv.xz})
      and rename \samp{DEW_aq.csv} to \samp{DEW.csv}.

      \item Remove old data file suffixes in references (e.g. [S92] was used to
      indicate that the data first appeared in \acronym{SUPCRT92}).

      \item Add GHS for almandine, dickite, fluorphlogopite, glaucophane,
      grunerite, larnite, pyrope (\samp{SUPCRT92.csv}) and bromellite
      (\samp{inorganic_cr.csv}) from
      \href{https://pubs.er.usgs.gov/publication/b2131}{Robie and Hemingway, 1995}.

    }
  }

  \subsection{DEMOS AND VIGNETTES}{
    \itemize{

      \item Add \samp{demo/comproportionation.R}: Gibbs energy of sulfur
      comproportionation, after
      \href{https://doi.org/10.1111/1462-2920.14982}{Amend et al., 2020}.

      \item Convert \samp{equilibrium.Rnw} to \samp{equilibrium.Rmd} and revise
      and simplify content.

      \item Revise \samp{OBIGT.Rmd} to reduce the size of the HTML file and
      make deep linking to individual sections work.

      \item Remove vignette \samp{hotspring.Rnw}, demo \samp{bison.R}, and
      associated data files \samp{DS11.csv} and \samp{DS13.csv}. The plots are
      now in the \samp{bison.Rmd} vignette in the \pkg{JMDplots} package
      (\url{https://github.com/jedick/JMDplots}).

    }
  }

  \subsection{CHANGES TO FUNCTIONS}{
    \itemize{

      \item \code{subcrt()}: change \strong{action.unbalanced} argument to
      \strong{autobalance}, which now provides the ability to prevent
      autobalancing.

      \item Changing the water model with \code{water()} now updates the
      literature references in \code{thermo()$OBIGT}.

      \item \code{thermo.refs()} now shows CHNOSZ version and date.

      \item \code{subcrt()} uses degree symbol in messages (\strong{°C}).

      \item Change \code{thermo$...} to \code{thermo()$...} in messages and
      comments.

    }
  }

  \subsection{OTHER CHANGES}{
    \itemize{

      \item Fix \samp{palply.Rd} for new warning about \dQuote{Non-file
        package-anchored link(s) in documentation object} in \command{R CMD
        check}.

      \item Convert this NEWS file to Rd format.

    }
  }

  \subsection{TODO}{
    \itemize{

      \item TODO: fix error in expr.species("AsO4---")

      \item TODO: make a site for community uploaded data

      \item TODO: refresh extdata/supcrt/compare.R

      \item TODO: get arsenopyrite from Perfetti et al., 2008

      \item TODO: make separate directories for default and optional data files

      \item TODO: demo/saturation.R: put Al+3 first in basis definition (it's
        confusing to have H2O first when the reactions are balanced on Al)

      \item TODO: diagram(): make error or warning for NA values of affinity

      \item TODO: use https in refs.csv

      \item TODO: move protein.equil() example (previously in equilibrium.Rnw) to JMDplots.

      \item TODO: change add.obigt --> add.OBIGT on website (Fe-O-H-S example) when package is submitted to CRAN.

      \item TODO: move extdata/bison to JMDplots package.

      \item TODO: remove extdata/refseq (now in JMDplots).

    }
  }

}

\section{Changes in CHNOSZ version 1.3.6 (2020-03-16)}{
  \itemize{

    \item Fix for R-devel (4.0.0): test-add.protein.R had a test that failed
      with new stringsAsFactors = FALSE default for read.table().

  }
}

\section{Changes in CHNOSZ version 1.3.5 (2020-02-16)}{
  \itemize{

    \item Fix for R-devel (4.0.0): remove documentation for unused argument
      '...' in add.protein.Rd.

  }
}

\section{Changes in CHNOSZ version 1.3.4 (2019-12-15)}{

  \subsection{THERMODYNAMIC DATA}{
    \itemize{
      \item Change data for dawsonite, scheelite, and ferberite from calories to
        Joules (as originally reported).

      \item Use correct Cp coefficients for ferberite (2nd term in equation of
        Wood and Samson, 2000 is off by a factor of 10). Thanks to Xiangchong
        Liu and David Polya for the bug report and assistance.

      \item Move more carbonates from SUPCRT92 back into the default database:
        artinite, azurite.
    }
  }

  \subsection{DIAGRAMS}{
    \itemize{
      \item Improve the handling of solids in equilibrate(). They are now excluded
        from the equilibrium calculation, but their stability fields are
        calculated using the maximum affinity method (with diagram()). Where
        solids are stable, logarithms of activities of aqueous species are set
        to -999, and vice versa. Thanks to Feng Lai for the request and test
        case.
    }
  }

  \subsection{PACKAGE CLEANUP}{
    \itemize{
      \item Remove files Sce.csv.xz and yeastgfp.csv.xz and functions yeastgfp()
        and yeast.aa(). These have been moved to the JMDplots package
        (https://github.com/jedick/JMDplots).

      \item Remove protein_refseq.csv.xz (based on RefSeq release 61). An updated
        version, based on RefSeq release 95, is available in JMDplots.
    }
  }

  \subsection{BUG FIXES}{
    \itemize{
      \item describe.reaction(): revert the Unicode double arrow. Unlike that, a
        simple equals sign is displayed correctly on all systems.

      \item Fix bug in unexported function obigt2eos(): lambda (exponent in heat
        capacity equation) was incorrectly going through a units conversion.

      \item subcrt(): properties of minerals are now output at the listed upper
        T limit (or transition temperature, for e.g. cr1 -> cr2). Previously,
        properties were set to NA at, and not only above, the T limit. Thanks
        to Evgeniy Bastrakov for the bug report.

      \item In mosaic(), set the activities of only aqueous basis species to the
        total activity, which is taken from the incoming basis() definition.
        This fixes a bug where activities of minerals (particularly, native
        sulfur) were unexpectedly changed; their logarithms are now set to 0.
        Thanks to Evgeniy Bastrakov for the bug report.

      \item Fix a bug where berman() returned NA values for K-feldspar below
        298.15 K, due to improperly initialized values for the disorder
        properties. Thanks to Kaustubh Hakim for the bug report.

      \item Use inherits() instead of class() for checking try() errors (don't
        use conditions of length > 1 in if() statements -- fix for R-devel).

      \item Fix a bug in affinity() with calculations in IS-T (but not T-IS)
        space. This involves a workaround for the limitation in subcrt() that
        only IS (if present) can be the 'grid' variable.
    }
  }

}

\section{Changes in CHNOSZ version 1.3.3 (2019-08-02)}{

  \subsection{NEW FEATURE: ENERGY UNITS}{
    \itemize{
      \item The csv data files in OBIGT now have a column named "E_units", which
        is used to specify the energy units for each species (cal or J).

      \item The messages produced by info() (and checkGHS() and checkEOS()) use
        cal or J depending on the units in the database. 

      \item All properties returned by info() reflect the units in the database
        for each species (no conversion is done).

      \item All properties returned by subcrt() continue to follow the user's
        units setting (from E.units(), default cal). The messages from
        subcrt() now include the units used for the output.

      \item For backward compatibility, add.obigt() allows reading data from a
        user-supplied file with no "E_units" column. In this case, the units
        are set to cal.
    }
  }

  \subsection{NEW FEATURE: SOLUBILITY UNITS}{
    \itemize{
      \item convert() now handles the output from solubility(), for converting
        calculated concentrations to ppm or ppb.

      \item diagram() has a new argument 'levels' for specifying the levels at
        which to draw solubility contours.

      \item demo/contour.R has been updated to show the 1, 10, 100, 1000 ppb
        contours for Au solubility.

      \item Add demo/sphalerite.R for calculations of solubility of sphalerite,
        converted to log ppm (after Akinfiev and Tagirov, 2014).
    }
  }

  \subsection{NEW FEATURES: DIAGRAMS}{
    \itemize{
      \item diagram(): add 'lty.cr' argument to specify line type for boundaries
        between minerals (two crystalline phases). Setting lty.cr = 0
        suppresses the lines, which is useful for plotting only the activities
        of aqueous species (can be used to make solubility diagrams for some
        systems).

      \item Also add 'lty.aq' to specify the line type between aqueous species.

      \item The 'names' argument of diagram() recognizes more values, making it
        easier to use in scripts. A value of NULL or TRUE plots the names, and
        NA or FALSE hides the names.

      \item Add functions lNaCl(), lS(), lT(), lP(), lTP(), and lex() for making
        legend text to describe system conditions.
    }
  }

  \subsection{BUG FIXES}{
    \itemize{
      \item In mosaic(), constant activities were used for each changing basis
        species, leading to an artificial inflation of the affinities of
        species around the transitions where the basis species become equally
        abundant. Now, affinities of formation of species are calculated using
        activities of basis species that are adjusted (i.e. lowered) according
        to their relative abundances.

      \item The 'mixing' argument of mosaic() has been removed; it is no longer
        needed with the improved handling of activities of basis species.

      \item Using the "Alberty" method, nonideal() calculated incorrect values of
        adjusted H, S, and Cp at given ionic strength. The correct equations
        are now used (taken from Alberty, 2001), and checks of the self-
        consistency of the calculated values and comparisons with the
        Debye-Hückel limiting slopes for G, H, and Cp have been added to
        test-nonideal.R.

      \item For improved calculation of the Debye-Hückel limiting slopes, the
        equation for the Debye-Hückel constant was changed to use the equation
        given by Clarke and Glew, 1980 (valid to 150 °C).
    }
  }

  \subsection{IMPROVEMENTS TO MOSAIC}{
    \itemize{
      \item In mosaic(), groups of changing basis species are equilibrated for a
        total activity given by the activity of the starting basis species in
        the incoming basis() definition. This provides a means to specify the
        the amount of a conserved element in the basis species.

      \item equilibrate(): make it possible to process the output of mosaic().
        This feature combines the equilibrium activities of the formed species
        with those of the changing basis species to make an object that can be
        plotted with diagram(). This can be useful for conserving two elements
        in a system: one in the basis species, the other in the formed
        species. Thanks to Kirt Robinson for the feature request and test
        system for acetamide, which has been added to ?equilibrate.

      \item Add tests to test-mosaic.R to check that activitives produced by
        mosaic() - equilibrate() and mosaic() - solubility() really are
        equilibrium activities, i.e. that the affinities of reactions between
        species are zero everywhere, in particular near the transitions of
        basis species.

      \item In mosaic(), include the relative abundances in each group of basis
        species in the output ('E.bases').
    }
  }

  \subsection{THERMODYNAMIC DATA}{
    \itemize{
      \item Correct the value of a1 for acetamide (Shock, 1993). Thanks to Kirt
        Robinson.

      \item Add properties of aqueous species (dimethylamine, trimethylamine,
        resorcinol, phloroglucinol, cyclohexane carboxylate, and cyclohexane
        carboxylic acid) from LaRowe and Amend, 2019 and HKF parameters for
        dimethylammonium chloride from Shvedov and Tremaine, 1997.

      \item Add data for zoisite, clinozoisite, and epidote from Gottschalk, 2004.
        Thanks to Evgeniy Bastrakov for contributing the data file.

      \item Move some carbonates from SUPCRT92 back into default database:
        huntite, hydromagnesite, malachite, nesquehonite.

      \item Add trisulfur radical ion (S3-) from Pokrovski and Dubessy, 2015.

      \item Correct entries for aqueous phenanthrene and methylphenanthrene
        isomers to use Joules (Dick et al., 2013).
    }
  }

  \subsection{DOCUMENTATION}{
    \itemize{
      \item demo/aluminum.R: add calculations using Si(OH)4 from the Akinfiev-
        Diamond model (SiO2 in these reactions is replaced by Si(OH)4 - 2 H2O).

      \item Modify demo/saturation.R to use the argument recall feature of
        affinity().

      \item In anintro.Rmd, include an example of using mod.obigt() to add a
        solid (magnesiochromite, with data from Klemme et al., 2000). Thanks
        to Alan Boudreau for the suggestion.

      \item In anintro.Rmd, add description of add.obigt("AkDi") and change
        mod.obigt() example for aqueous species to use CoCl4-2 (Liu et al.,
        2011).
    }
  }

  \subsection{OTHER CHANGES}{
    \itemize{
      \item Add function moles(), to calculate the total number of moles of
        elements in the output from equilibrate().

      \item subcrt(): stop with an error if length(logact) is not the same as the
        number of species.

      \item Add 'transect' argument to affinity(), to allow forcing transect mode
        for any number of values of the variables.

      \item solubility() has a new argument 'codeanal' that when set to TRUE
        indicates the printing of additional diagnostic messages.

      \item check.obigt(): add 'E_units' column and fix a bug, also involving
        makeup() and entropy(), where many species with differences between
        database and computed values of Gibbs energy were not listed.
    }
  }

}

\section{Changes in CHNOSZ version 1.3.2 (2019-04-20)}{

  \subsection{NEW FEATURES}{
    \itemize{
      \item Add thermo/stoich.csv.xz (loaded as thermo()$stoich), containing a
        precalculated stochiometric matrix for the default database, to speed
        up retrieve().

      \item retrieve() now updates the stoichiometric matrix when the database
        changes, instead of performing a full recalculation.

      \item Add 'ligands' argument to retrieve(), for getting metal-bearing
        species with a range of possible elements in the ligands.

      \item diagram(): Starting with R version 3.6.0, 'fill' can be the name of
        any available HCL-based color palette (see hcl.pals() for a list).

      \item Add demo/contour.R for gold solubility contours on logfO2-pH diagram.

      \item In the output of affinity(), the list of 'vals' (the values of the
        input variables) now has names corresponding to the variables.

      \item In diagram(), a background grid corresponding to major and/or minor
        tick marks can be added using the 'grid' argument. Or, use
        thermo.axis() after diagram() to overlay a grid.

      \item Add demo/arsenic.R (Eh-pH diagram after Lu and Zhu, 2011).

      \item Rename demo/go-IU.R to aluminum.R, and add plots from Tutolo et al.,
        2014.

      \item Add argument 'plot.box' (default TRUE) to thermo.plot.new() and
        'plot.line' (default FALSE) to thermo.axis() to prevent unnecessary
        plotting of axis lines. thermo.axis() also avoids plotting tick
        marks at the plot limits.
    }
  }

  \subsection{THERMODYNAMIC DATA}{
    \itemize{
      \item OBIGT: use dawsonite Cp coefficients from Tutolo et al., 2014.

      \item OBIGT: Move H4SiO4 (Stefánsson, 2001) to AS04.csv (optional data).

      \item Berman: rename .csv files to include whole year. berman()
        now uses the years in the filenames to assemble the data in reverse
        chronological order, and keeps only the most recent data for any
        duplicated minerals.

      \item Berman: Add data for sudoite, daphnite, Mg-amesite, and Fe-amesite
        (Vidal et al., 1992, 2001, 2005).

      \item Berman: make "albite" identical to "albite,low", not "albite,high".

      \item Move transuranic actinide complexes to optional data (SLOP98), as
        entropies of the elements at 298.15 K are not available to check the
        self-consistency of the GHS values for the complexes.

      \item obigt.Rmd: Use Zenodo DOI for citations of slop98.dat, slop07.dat, and
        slop16.dat.
    }
  }

}

\section{Changes in CHNOSZ version 1.3.1 (2019-03-02)}{
  \itemize{

    \item Fix unexported function obigt2eos(), which was affected by "Assignment
      to empty subset of data.frame fails check of length of value"
      (https://bugs.r-project.org/bugzilla/show_bug.cgi?id=17483). The bug
      caused errors in R-release (3.5.2) but not R-devel, when using info()
      for species with the AkDi model.

    \item Fixes to maintain compatibility with R version 3.1.0: don't use
      lengths(), and provide explicit colClasses to read.csv() to avoid
      automatic assignment of character for long numeric values.

  }
}

\section{Changes in CHNOSZ version 1.3.0 (2019-02-26)}{

  \subsection{CRAN COMPLIANCE}{
    \itemize{
      \item Because data(thermo) and data(OBIGT) did not actually create said
        objects in the user's workspace (the expected behavior for data()),
        the package's data are now loaded by two new functions, reset() and
        obigt(). reset() is automatically run upon loading the package, so
        that existing scripts beginning with data(thermo) still work (this
        command now has no effect other than producing a warning).
    }
  }

  \subsection{NEW FEATURES}{
    \itemize{
      \item Add thermo() as a convenience function to access or modify the
        package's data, especially various computational options.

      \item Add retrieve() to retrieve all the species having given elements, or
        all species in a given chemical system. Thanks to Evgeniy Bastrakov
        for the suggestion.

      \item Add AkDi() to calculate thermodynamic properties of aqueous
        nonelectrolytes using the Akinfiev-Diamond model. Thanks to Evgeniy
        Bastrakov for guidance.
    }
  }

  \subsection{THERMODYNAMIC DATA}{
    \itemize{
      \item Revert to using SiO2(aq) from SUPCRT92 (i.e. Shock et al., 1989) in
        the default database.

      \item Move SiO2(aq) from Apps and Spycher, 2004 and recalculated HSiO3- to
        new optional data file, OBIGT/AS04.csv.

      \item Move H2AsO3- from OBIGT/SLOP98.csv to OBIGT/SUPCRT92.csv.

      \item Add OBIGT/AkDi.csv as optional data file for aqueous species from
        Akinfiev and Diamond, 2003 and Akinfiev and Plyasunov, 2014.

      \item Add gaseous HF and HCl to as requirements for the Akinfiev-Diamond
        model for dissolved gases.

      \item Move Au+3, AuCl3-2, and AuCl4-3 back into default database.

      \item Add/replace data for Ag+ and Cu+ and their complexes from Akinfiev
        and Zotov, 2001 and 2010.

      \item Update HCl(aq) with data from Tagirov et al., 1997.

      \item Update Pd+2 and complexes with data from Tagirov et al., 2013.

      \item Update Pt+2 and complexes with data from Tagirov et al., 2015.

      \item Update Zn+2 and complexes with data from Akinfiev and Tagirov, 2015.
    }
  }

  \subsection{DOCUMENTATION}{
    \itemize{
      \item In demo/NaCl.R, indicate region not considered by Shock et al., 1992
        in developing the "g function" applicable to electrolytes in the
        revised HKF model, and note presence of resulting discontinuities (see
        man/examples.Rd).

      \item Add demo/AkDi.R showing calculations of Henry's constant of dissolved
        gases, after Figure 1 of Akinfiev and Diamond, 2003.
    }
  }

  \subsection{OTHER CHANGES}{
    \itemize{
      \item In describe.reaction(), change equals sign to reaction double arrow.

      \item water() and water.SUPCRT92(): add 'P1' argument to choose whether to
        output 1 bar for Psat at temperatures less than 100 degrees C
        (default is TRUE).

      \item diagram(): 'xlim' and 'ylim' arguments now apply to 2-variable
        diagrams. Thanks to Evgeniy Bastrakov for the suggestion.
    }
  }

}

\section{Changes in CHNOSZ version 1.2.0 (2019-02-09)}{

  \subsection{BUG FIXES}{
    \itemize{
      \item Fix bug in nonideal() where "Zn" in formula was identified as charge.
        Thanks to Feng Lai for the bug report.

      \item Fix a bug where subcrt()$reaction$coeffs was incorrect for reactions
        involving minerals with phase transitions. Also ensure that the output
        reaction stoichiometry is correct for duplicated species in reactions.
        Thanks to Grayson Boyer for the bug report.

      \item For species in the revised HKF model, subcrt() now sets properties to
        NA where the density of H2O is less than 0.35 g/cm3, avoiding the
        output of bogus values in this region. Thanks to Evgeniy Bastrakov.

      \item For systems where no basis species is present in all formation
        reactions, and the user hasn't provided balance coefficients, stop
        with an error instead of setting the balance cofficients to 1.
        Thanks to Shuang Kong for an example calculation and Tucker Ely for
        a previous suggestion to produce an error here. The affected code is
        in balance(), an unexported function used in equilibrate() and
        diagram() (and now also solubility()).
    }
  }

  \subsection{NEW FEATURE: SOLUBILITY CALCULATIONS}{
    \itemize{
      \item Add solubility(). Run this after affinity() to calculate the
        solubility of a solid or gas defined as the conserved basis species,
        which is involved in the formation of one or more dissolved species.

      \item Features include automatic detection of dissociation reactions and
        finding the final ionic strength for dissolution of a mineral into
        pure water (find.IS argument).

      \item find.IS depends on the new argument recall feature of affinity() (or
        mosaic() if that is used instead). This allows a calculation to be
        re-run with the same settings except for particular additions or
        modifications, in this case modified values of ionic strength.

      \item Revise demo/solubility.R to show solubility calculations for CO2(gas)
        and calcite as a function of T and pH.
    }
  }

  \subsection{NEW FEATURE: EXPANDED ACTIVITY COEFFICIENT CALCULATIONS}{
    \itemize{
      \item Reorganize and expand options for activity coefficient calculations
        (set in thermo$opt$nonideal: Bdot, Bdot0, bgamma, bgamma0, or Alberty).
        The previous default, which corresponds to 'bgamma' (T- and
        P-dependent extended term parameter with single ion-size parameter),
        has been replaced by 'Bdot' (T-dependent extended term parameter and
        species-dependent ion-size parameter; see below).

      \item nonideal() with the 'Bdot' or 'Bdot0' equation uses specific
        ion-size parameters for different ions, in accord with the HCh package
        (Shvarov and Bastrakov, 1999). Parameters are from Table 2.7 of
        Garrels and Christ, 1965.

      \item nonideal() now calulates activity coefficients of neutral species,
        using the Setchénow equation. Whether the extended-term parameter in
        this equation is taken to be zero or is taken from the value for
        charged species (see above) is controlled by setting
        'thermo$opt$Setchenow' to bgamma0 (default) or bgamma.

      \item Add NaCl(), implementing a first-order calculation of the speciation
        of NaCl in water, taking account of activity coefficients and the
        reaction Na+ + Cl- = NaCl(aq).
    }
  }

  \subsection{OTHER NEW FEATURES}{
    \itemize{
      \item Add dumpdata() for returning/writing all packaged thermodynamic data
        (including default database and optional data files). The file is
        also available on the website (chnosz.net/download/alldata.csv).

      \item mosaic() has been rewritten to handle more than two changing groups
        of basis species.

      \item mosaic() gets a new argument 'mixing' (default TRUE), indicating
        whether the Gibbs energy of ideal mixing should be included in the
        calculations with blend = TRUE. As a test of this change, an example
        using mosaic() to calculate the pH-dependent solubility of calcite
        has been added to solubility.Rd.
    }
  }

  \subsection{DOCUMENTATION}{
    \itemize{
      \item Add demo/gold.R for calculations of Au solubility in hydrothermal
        chloride and sulfide solutions (based on diagrams from Akinfiev and
        Zotov, 2001, Stefánsson and Seward, 2004, and Williams-Jones et al.,
        2009). This depends on the revised nonideal() and new NaCl() functions
        described above.

      \item anintro.Rmd: add section on solubility calculations. Add cuprite to
        mosaic diagram example, and note about implications of changing
        balance coefficients to 1. Add 'mklinks.sh' script for adding links
        to CHNOSZ and R documentation.

      \item Add demo/bison.R (average oxidation state of carbon of metagenome-
        derived proteins in different microbial phyla at Bison Pool)

      \item demo/saturation.R replaces the previous activity_ratios.R, and
        includes not only activity ratios on the axes but also saturation
        limits for calcite, dolomite, magnesite, and brucite in the
        H2O-CO2-CaO-MgO-SiO2 system.

      \item Add demo/glycinate.R showing logK of complexation of glycinate with
        divalent and monovalent metals.

      \item Remove vignette wjd.Rnw (Winding journey down in Gibbs energy).
    }
  }

  \subsection{THERMODYNAMIC DATA}{
    \itemize{
      \item The Berman data (Berman, 1988 and later additions) have replaced the
        SUPCRT92 data (based on Helgeson et al., 1978) for most minerals in
        the default database (i.e. the one loaded by data(thermo)). Only
        native elements, sulfides, halides, sulfates, and selected carbonates
        and oxides from SUPCRT92 that do not duplicate any in the Berman
        dataset are loaded by default.

      \item The superseded SUPCRT92 minerals (primarily the silicates and
        Al-bearing minerals, as well as calcite, dolomite, hematite, and
        magnetite) are in the file extdata/OBIGT/SUPCRT92.csv, which can be
        loaded with add.obigt("SUPCRT92"). This is used in some examples for
        comparing the datasets (see ?berman and demo("go-IU")), but should not
        be needed by most users.

      \item In thermo$obigt, the minerals using the Berman equations are
        identified as those having NA values for all thermodynamic properties
        (columns 8-20). These minerals have the state of "cr" (not "cr_Berman"
        as before).

      \item The heat capacity equation in berman() has been expanded to
        Cp = k0 + k1T-1/2 + k2T-2 + k3T-3 + k4T-1 + k5T + k6T2 (as used in
        winTWQ).

      \item Multipliers on values in the Berman data files (extdata/Berman/*.csv)
        have been removed, except for 10^5, 10^5, 10^5, and 10^8 on v1 to v4,
        to be consistent with the winTWQ data format.

      \item berman() now reads data from a user-supplied file specified in
        thermo$opt$Berman.

      \item Add mineral data from Berman and Aranovich, 1996 to
        extdata/Berman/testing and an example of loading these data via
        add.obigt() and thermo$opt$Berman (see ?berman).

      \item Correct swapped v1 <-> v3 and v2 <-> v4 in extdata/Berman/DS10.csv.

      \item To reduce confusion that "d5" might be a coefficient in the same
        equation as d0 .. d4, rename "d5" (constant that scales disordering
        enthalpy to volume) to "Vad" (notation used in Theriak-Domino manual).

      \item Change abbreviation of grossular to Grs.

      \item Make names for dicarboxylic acid anions from Shock, 1995 (in
        organic_aq.csv) consistent with those used in later papers, e.g.
        H-succinate and succinate-2 for the -1 and -2 charged ions. Remove
        duplicated entries for H-succinate and succinate-2 from biotic_aq.csv.

      \item Update G, H, and S of celestite using properties of the dissolution
        reaction given by Reardon and Armstrong, 1987.

      \item Add scheelite (CaWO4) and ferberite (FeWO4) using data adopted by
        Wood and Samson, 2000.

      \item Update aqueous Au species with data from Akinfiev and Zotov, 2001
        and 2010, and Pokrovski et al., 2014.

      \item Move SUPCRTBL updates (As and Al minerals and aqueous species and
        SiO2(aq)) into default database.

      \item Values for SiO2(aq) are taken from Apps and Spycher, 2004. Thanks
        to John Apps for providing the document.

      \item Recalculate GHS for HSiO3- by adding difference from SiO2(aq)
        (Sverjensky et al., 1997) to updated values for SiO2(aq).

      \item Add test-recalculate.R to check that recalculated values are correctly
        entered in OBIGT.

      \item Move superseded data for aqueous Au, As, and Al species, SiO2 and
        HSiO3- to OBIGT/SLOP98.csv.

      \item Add or update glycine, diglycine, and triglycine (zwitterions and
        ions), and diketopiperazine, [Gly] and [UPBB] groups from Kitadai,
        2014.

      \item Add data for metal-glycinate complexes from Azadi et al., 2019,
        superseding Shock and Koretsky, 1995.

      \item For reproducing previous calculations, superseded data for [Gly] and
        [UPBB], as well as [Met] (earlier superseded by LaRowe and Dick,
        2012), and metal-glycinate complexes have been moved to
        OBIGT/OldAA.csv.
    }
  }

  \subsection{DIAGRAMS}{
    \itemize{
      \item Lines in 1-D diagram()s can optionally be drawn as splines using the
        method for splinefun() given in the 'spline.method' argument (the
        default of NULL means no splines).

      \item Add 'srt' argument to diagram() (rotation of line labels).

      \item Export thermo.axis(), as it is useful for adding major and minor tick
        marks after (above) other plot elements such as legends.

      \item diagram(): rename 'what' argument to 'type'.

      \item diagram(): add new type of diagram, 'saturation', which is used to
        plot saturation lines for minerals (where their affinity equals
        zero).

      \item thermo.plot.new() and, by extension, diagram(), now "listen" to the
        setting of par("mar"). However, if a new plot is being made, or the
        device has the default setting of mar from par(), mar is (re)set to
        the default value in thermo.plot.new().

      \item diagram(): numeric 'names' now indicates subset of species names to
        include as labels on the plot.

      \item diagram(): new argument 'family' to select the font family for species
        names (labels on the plot). New arguments 'bold' and 'italic' apply to
        formatted names (e.g. chemical formulas), and 'font' applies to those
        without formatting.

      \item 'lines' is added to the return list of diagram(), giving the
        coordinates of lines (field boundaries) on 2-D diagrams (these are
        taken from the output of contourLines()).
    }
  }

  \subsection{COMPUTATIONAL OPTIONS}{
    \itemize{
      \item Add 'exceed.rhomin' argument to subcrt() and affinity() to enable
        output of properties for species in the revised HKF model below 0.35
        g/cm3.

      \item In equilibrate(), accept a length > 1 'normalize' argument to
        normalize the chemical formulas of only the selected species.

      \item Add thermo$opt$maxcores (default 2) to specify maximum number of
        cores for parallel calculations with palply().

      \item Keywords in basis(): Change 'CHNOPS+' to use O2 instead of e-, and add
        'CHNOPSe' and 'MgCHNOPSe' for sets of basis species that have e-.

      \item Add 'keep.duplicates' argument to thermo.refs(). Set it to TRUE to
        output a single primary reference for each species, keeping any
        duplicated references (but not including any secondary references in
        thermo$obigt$ref2). Thanks to Evgeniy Bastrakov for the suggestion.
    }
  }

  \subsection{USABILITY ENHANCEMENTS}{
    \itemize{
      \item To provide better diagnostics for other apps using CHNOSZ, warning
        messages produced by subcrt() are now available in the output of
        affinity(), under 'sout$warnings'.

      \item Change internal variable names in subcrt() for better readability
        (sinfo -> ispecies, inpho -> iphases, sinph -> phasespecies).

      \item info(): display mineral names instead of just "cr" in messages (e.g.
        for info("SiO2")).

      \item info(): species with no matches in default database are also
        searched in optional data files.

      \item Add C implementation of counting occurrences of all letters in a
        string (src/count_letters.c) to speed up operation of count.aa().

      \item read.fasta(): add support for file connections created using
        archive::archive_read (https://github.com/jimhester/archive).

      \item The arguments in expr.species() have been reorganized for more
        flexible and concise usage.

      \item Reduce message clutter for diagrams showing water stability limits.

      \item In subcrt(), improve error message for invalid property names.

      \item In affinity(), make invalid variable names an error.

      \item Add 'return.sout' argument to affinity(), to return just the values
        calculated with subcrt().
    }
  }

}

\section{Changes in CHNOSZ version 1.1.3 (2017-11-13)}{

  \subsection{MAJOR CHANGES}{
    \itemize{
      \item Add berman() function and extdata/Berman/*.csv files for calculating
        thermodynamic properties of minerals using equations of Berman, 1988.

      \item Calculations related to Berman's (1988) Figs. 1 and 2 for the lambda
        transition of quartz are available in the new demo lambda.R.

      \item Add functions implementing the Deep Earth Water (DEW) model
        (Sverjensky et al., 2014): water.DEW() and its supporting functions
        calculateDensity(), calculateGibbsofWater(), calculateEpsilon(),
        calculateQ().

      \item The computational setting for water (thermo$opt$water) can now be set
        using water("DEW"), water("IAPWS"), etc.

      \item Usage of the DEW model is shown in the new demo DEW.R. This demo also
        depends on the Berman equations (above) and, for the last diagram in
        the demo, the following *four* NEWS items:

      \item In equilibrate(), it is now possible to combine affinity calculations
        with variable activity of the balancing basis species (loga.balance).
        For example, in the last plot of the DEW demo, the calculation
        involves simultaneously varying temperature and logfO2 (given as
        arguments to affinity()) as well as total concentration of carbon
        (given by the loga.balance argument in equilibrate()).
    }
  }

  \subsection{COMPUTATIONAL IMPROVEMENTS}{
    \itemize{
      \item The 'alpha' argument of diagram() can be set to 'balance' to scale
        the values by the balancing component. This is useful for making
        "percent carbon" plots for systems where the species have different
        carbon numbers.

      \item nonideal() now has three methods for calculating activity
        coefficients: 'Helgeson' and 'Helgeson0', which utilize the
        Debye-Huckel equation with parameters as described by Helgeson, 1969,
        and 'Alberty' (the only method previously available). The 'Helgeson'
        method depends on the following *two* NEWS items (but the 'Helgeson0'
        method omits the B-dot term):

      \item All three water options (SUPCRT92, IAPWS95, DEW) can be used to
        calculate 'A_DH' and 'B_DH', i.e. A and B coefficients in the extended
        Debye-Huckel equation as described by Helgeson, 1969.

      \item Add Bdot() to calculate the "B-dot" (or b_gamma) extended term
        parameter for activity coefficients in NaCl solutions at high
        temperature and pressure (Helgeson et al., 1981), including
        very high pressures (> 10 kbar) based on extrapolations from Manning
        et al., 2013.

      \item For minerals with phase transitions (states 'cr2' 'cr3' etc.) in
        thermo$obigt (i.e. the Helgeson minerals), it is now possible to use
        the minerals in basis(), species(), affinity() with proper accounting
        for transition effects. Just set the state to 'cr' or omit it,
        e.g. species("quartz").

      \item Implement SUPCRT92's handling of variable volume for quartz and
        coesite. Calculations for other minerals still assume constant
        volume of each phase.

      \item Calculations of the g function are now enabled for DEW (with
        pressure derivative) and IAPWS-95 (no derivatives included).

      \item water.lines() now works for diagrams of Eh, pe, logfO2, logaO2,
        logfH2, or logaH2 vs pH, T, or P. It is possible to have T or P on
        either the x- or y-axis.

      \item Variables that are in the aguments to affinity() now override any
        buffers previously specified using basis().
    }
  }

  \subsection{DATABASE UPDATES}{
    \itemize{
      \item Add data(OBIGT) command to reset only the thermodynamic database
        (OBIGT) to its default. data(thermo) is still used to reset
        *everything* (the database and all computational and system settings).

      \item In add.obigt(), using 'file' without a file suffix now can be used to
        locate corresponding files (such as DEW_aq.csv or SUPCRTBL.csv) in the
        extdata/OBIGT package directory.

      \item add.obigt() gets new argument 'species' for selecting species to add.

      \item Add Berman_cr.csv with names and formulas for minerals whose
        thermodynamic properties are calculated using the Berman, 1988
        equations. The state "cr_Berman" is used to distinguish these
        minerals from the Helgeson data. The parameters themselves are
        stored in extdata/Berman/*.csv.

      \item Add DEW_aq.csv with aqueous species data from the DEW spreadsheet
        (May 2017 version). Species with data that do not differ from the
        default database in CHNOSZ are not included. Currently, detailed
        reference keys and descriptions (in thermo$refs) are given for
        selected species only. [Optional Data; not loaded by data(thermo)]

      \item Add SUPCRTBL.csv with mineral and aqueous species data compiled by
        Zimmer et al., 2016. Data were taken from original references or
        from spronsbl.dat (http://www.indiana.edu/~hydrogeo/supcrtbl.html).
        The new demo/go-IU.R uses these data. [Optional Data]

      \item Update HKF parameters of adenine based on experimental heat capacity
        and volume data (Lowe et al., 2017). Add demo/adenine.R comparing
        the updates with the previous data (LaRowe and Helgeson, 2006).

      \item Remove H2O(gas) from the database (taken from SUPCRT92), as it gives
        inconsistent results for metastable steam (Zavarin et al., 2016,
        LLNL-TR-701407, doi: 10.2172/1325873).

      \item Change the formula of the electron from "Z0-1" (parsed as "zero Z
        with minus one charge") to the equivalent but nicer-to-read "(Z-1)"
        (parsed as "minus one Z"; the parentheses cause "-1" to be interpreted
        as the coefficient on Z instead of the charge in the formula).

      \item Add data for hedenbergite, andradite, ferro-actinolite, grunerite, and
        ilvaite from Delgado Martin and Soler i Gil, 2010 (using the Berman
        equations).
    }
  }

  \subsection{OTHER CHANGES}{
    \itemize{
      \item Add test-logmolality.R to demonstrate transformation of variables from
        activity to molality in the main workflow; refer to this test and
        describe these transformations in anintro.Rmd; remove "experimental"
        labeling of activity coefficient calculations in CHNOSZ-package.Rd.

      \item Add 'use.molality' argument to axis.label() and supporting functions.
        So that the labels on diagrams are changed from activity (a) to
        molality (m), this argument is set to TRUE by diagram() when IS is
        supplied as an argument to affinity().

      \item Add demo TCA.R for standard Gibbs energies of steps of the
        citric acid cycle (Canovas and Shock, 2016).

      \item Components of subcrt() output indicating the stable polymorph of
        minerals with phase transitions are now named 'polymorph', not
        'state'.

      \item Add maxdiff() and expect_maxdiff() for calculating and testing the
        maximum absolute pairwise difference between two objects.
        expect_maxdiff() is used in tests/testthat/*.R as a more
        straightforward alternative for checking differences than using
        'tolerance' in expect_equal().

      \item Add 'tol' argument to equil.reaction() (convergence tolerance for
        uniroot()).

      \item Add thermo$opt$Berman to signal info(), subcrt(), etc. to prefer
        minerals in Berman dataset over Helgeson dataset (default: FALSE).
    }
  }

  \subsection{REFACTORING AND CLEANUP}{
    \itemize{
      \item In hkf() and cgl(), combine 'ghs' and 'eos' arguments into single
        argument named 'parameters'.

      \item Remove 'H2O.PrTr' and 'domega' arguments from hkf(). Properties of
        water at Pr and Tr and the logic for when to find derivatives of omega
        are now calculated within the function, thereby simplifying the
        function call from subcrt().

      \item Add 'H2O.props' argument to hkf(). This lists the properties of water
        to calculate (at P, T); these are now included in the results, in
        order to avoid duplicating the calculations in subcrt().

      \item In water.* functions, rename the "diel" variable to "epsilon".

      \item To save space, taxid_names.csv has been trimmed to hold only those
        taxids listed in extdata/bison/gi.taxid.txt.

      \item Remove Eco.csv.xz and rename more.aa() to yeast.aa().

      \item Remove transfer() and associated functions draw.transfer(), apc(), and
        feldspar().

      \item Remove read.expr() and extdata/abundance/ISR+08.csv (protein abundance
        in E. coli).

      \item Reduce the size of obigt.html by ~450 KB. This was done by 1) using
        rmarkdown with files in rmd/h/bootstrap-3.3.5/fonts/ replaced with
        empty files (291 KB), and 2) changing the theme from "cosmo" to
        "default" (151 KB).

      \item Remove extdata/OBIGT/CHNOSZ_aq.csv and move the entry for
        pseudo-H4SiO4 to inorganic_aq.csv.

      \item Reorganize concept index entries in Rd files: Main workflow, Extended
        workflow, Thermodynamic data, Thermodynamic calculations, Water
        properties, Protein properties, Other tools, Utility functions.

      \item Remove anim.TCA() and anim.carboxylase(); move the latter code to
        demo/carboxylase.R.
    }
  }

}

\section{Changes in CHNOSZ version 1.1.0 (2017-05-04)}{

  \subsection{DOCUMENTATION}{
    \itemize{
      \item Replace anintro.Rnw (Sweave) with rewritten anintro.Rmd (using knitr
        with Tufte style).

      \item Add obigt.Rmd (Thermodynamic data in CHNOSZ).

      \item Add eos-regress.Rmd (Regressing thermodynamic data), including
        section on Making a pseudospecies: H4SiO4.

      \item New demos: bugstab.R (potential diagrams for microbial proteins in
        colorectal cancer), activity_ratios.R (mineral stability diagrams
        with activity ratios on the axes), Shh.R (affinities of
        transcription factors relative to Sonic hedgehog).

      \item Converted demos: move examples from help pages to protein.equil.R,
        add.obigt.R, and affinity.R.

      \item Updated and removed demos: update buffer.R (logfH2 set by mineral
        buffers and aqueous species) and rename old buffer.R to protbuff.R
        (proteins as buffers); remove nucleobase.R.
    }
  }

  \subsection{NAMESPACE}{
    \itemize{
      \item The NAMESPACE now exports selected functions, instead of all of them.

      \item The following functions are no longer exported, nor documented in the
        Rd files: thermo.axis, optimal.index, extremes, caller.name, lsub,
        lsum, pprod, energy, energy.args, A.ionization, can.be.numeric, s2c,
        c2s, eos.args, TP.args, state.args, info.approx, info.character,
        info.text, info.numeric, count.charge, count.formulas, buffer,
        mod.basis, put.basis, preset.basis, preset.logact, species.basis,
        envert, outvert, which.balance, gfun, IAPWS95.residual,
        IAPWS95.idealgas, obigt2eos, group.formulas, get.formula, get.objfun,
        balance.
    }
  }

  \subsection{CHANGES TO diagram()}{
    \itemize{
      \item For the `groups` argument, activities are multiplied by the balance
        coefficients before the summation.

      \item Draw smoother boundary lines (with contour(); `lty` is usable), at the
        expense of not being precisely aligned with the grid used for plotting
        colors. For the old behavior, use `dotted=0`.

      \item Default of `legend.x` changed to NA (instead of making a legend, put
        labels near the lines).

      \item Field labels: apply formatting (subscripts and superscripts) to names
        that parse as chemical formulas; use `format.names=FALSE` to disable.

      \item Chemical formulas that are part of the axis labels (log activity or
        fugacity) now have formatting applied.

      \item Add arguments `adj` and `dy` for x-alignment and y-offset of line
        labels.

      \item Add arguments `fill.NA` (color of empty areas) and `limit.water`
        (assign NA to areas beyond water stability limits on Eh-pH and pe-pH
        diagrams).
    }
  }

  \subsection{NEW FEATURES}{
    \itemize{
      \item Add ZC.col() for generating a red-grey-blue color scale from
        numeric values.

      \item Add "QEC" as a keyword for preset species in basis() (glutamine,
        glutamic acid, cysteine, H2O, O2).

      \item More flexible parsing of chemical formulas for ZC() and other
        functions; e.g. `ZC(colMeans(protein.formula(1:4)))` now works.

      \item read.fasta() extracts information from UniProt FASTA headers
        (accession, name, organism) into columns of the output data frame
        (protein, abbrv, organism).

      \item Add `...` argument to thermo.plot.new() (additional arguments for
        par()).

      \item With new default arguments, thermo.axis() only plots the axis tick
        marks and lines; this is used in diagram() to redraw the axes on
        filled diagrams.

      \item seq2aa() removes newlines and whitespace before counting the
        letters in the sequence.

      \item read.fasta(): change argument `i` to `iseq`; this is used to select
        particular sequences to read from the file.

      \item Add syslab() and ratlab() to make formatted expressions for
        sets of thermodynamic components and activity ratios.

      \item Add GHS_Tr() to calculate G, H, and S at Tr (298.15 K) of phases
        stable at high temperature, using provided enthalpies of transition
        and database values of transition temperatures and heat capacity
        parameters, and GHS at Tr of the low-temperature phase.

      \item Add extdata/supcrt/ directory with files read.supcrt.R, compare.R,
        newnames.csv; used to compare data in SUPCRT and slop files
        (http://geopig.asu.edu/?q=tools) to data in thermo$obigt.
    }
  }

  \subsection{NEW AND OLD DATA}{
    \itemize{
      \item Split OBIGT.csv and OBIGT-2.csv into multiple files according to type
        of substance (inorganic, organic, or biotic). All files are now under
        extdata/OBIGT/. Conflicting data from OBIGT-2.csv were removed (see
        below).

      \item Citations in OBIGT now have a suffix (e.g. [S92] or [S98]) if the data
        were first present in sprons92.dat from SUPCRT, slop98.dat, or later
        versions of the slop data files (http://geopig.asu.edu/?q=tools).

      \item biotic_aq.csv: Add data for citric acid cycle metabolites from
        Canovas and Shock, 2016.

      \item CHNOSZ_aq.csv: Add H4SiO4 pseudospecies (pseudo-H4SiO4); see
        calculations in eos-regress.Rmd.

      \item organic_aq.csv: Add data for aqueous phenanthrene and methylphenanthrene
        from Dick et al., 2013.

      \item Remove superseded data previously in OBIGT-2.csv: citric acid cycle
        metabolites from Dalla-Betta and Schulte, 2009.

      \item Remove superseded data previously in OBIGT-2.csv: methionine and [Met]
        from Dick et al., 2006. Examples that require the old group additivity
        parameters can modify the database using:
        `mod.obigt("[Met]", G=-35245, H=-59310)`.

      \item Move two datasets that conflict with the primary database, but are
        used in tests, examples, or vignettes, from OBIGT-2.csv to
        extdata/thermo (Ste01.csv and BZA10.csv).

      \item Remove other conflicting datasets (i.e. that provide different data
        for a species already in the database) previously listed in
        OBIGT-2.csv (21 sources).

      \item Replace GHS for iron (c2, c3) calculated in CHNOSZ with values derived
        from sprons92.dat.

      \item Replace data for sulfur (cr) from Gurrieri, 1996, personal
        communication with that "as used by" McCollom and Shock, 1997 (from
        slop98.dat).

      \item Update selected nonelectrolytes with data from Plyasunov and Shock,
        2001, as listed in slop07.dat (http://geopig.asu.edu/?q=tools).

      \item extdata: Add fasta/rubisco.fasta and cpetc/*.csv files SS98_Fig5a,
        SS98_Fig5b, PM90 (all used in anintro.Rmd).
    }
  }

  \subsection{BUG FIXES}{
    \itemize{
      \item Correct charge (-2) of NAD(red)-2 in biotic_aq.csv. Thanks to Peter
        Canovas.

      \item Correct charge (-2) of MgATP-2.

      \item subcrt() returns `loggam` using the common logarithm; add
        test-nonideal.R to check for consistency between loggam and logK
        values returned by subcrt(). This also fixes the issue previously
        noted for the first example in ?nonideal. Thanks to David T. Wang for
        the bug report and test.

      \item NaN values from equilibrate() are now preserved by diagram(),
        producing empty (NA) fields rather than being mistakenly labeled with
        the first species. Thanks to Grayson Boyer for the bug report.
    }
  }

  \subsection{OTHER CHANGES}{
    \itemize{
      \item Remove msgout(), and replace previous calls to that function with
        message() from base R. As a result, the messages don't appear in
        Sweave vignettes, but can be turned on or off in knitr vignettes.

      \item Remove HTCC1062.faa.xz.

      \item Remove the bold formatting for "A" in expr.property(); "A" is
        italicized like most other uppercase letters.

      \item Rename browse.refs() to thermo.refs(); remove URL browsing (except for
        summary table).

      \item New function pinfo() merges functionality of old iprotein() and
        ip2aa(). Add `regexp` argument to control whether matches are made
        using a regular expression.

      \item Rename aa2eos() to protein.obigt().

      \item Remove anim.plasma() and AA03.csv.

      \item Remove stress() and stress.csv; move data from Tai et al., 2005 (used
        in an example in ?read.expr) to TBD+05.csv.

      \item Remove read.aa() - replaced by read.csv() with as.is=TRUE.

      \item Remove grep.file().

      \item Remove sideeffects.Rd (functions with side effects are now identified
        with red color in anintro.Rmd).
    }
  }

}

\section{Changes in CHNOSZ version 1.0.8 (2016-05-28)}{
  \itemize{

    \item Add gypsum to OBIGT.csv.

    \item Add MgSO4 and hydrated Mg sulfate minerals to OBIGT-2.csv (kieserite,
      starkeyite, hexahydrite and epsomite from Grevel and Majzlan, 2009).

    \item Add demo dehydration.R. This demo requires the RSVGTipsDevice package,
      which is currently unavailable on Windows.

    \item Add ... argument to EOSvar() and EOSregress(). Allows specifying new
      variables that are written as a function of arbitrary properties
      (regression variables still must be a function of at least T and P).

    \item EOSvar() now returns the negative of the 'QBorn' variable so that
      values of omega obtained from regressions of volume data have the
      correct sign in the HKF equations.

    \item Use given/family naming scheme in inst/CITATION.

    \item makeup() now works even if 'thermo' object is not available (useful
      for other packages importing this function). Thanks to Nir Shahaf for
      the suggestion.

    \item So that byte-compiled package code works, subcrt() no longer uses
      caller.name() for flow control. In its place, a new argument
      'check.Ttr' has been added, used to bypass the checks of transition
      temperature in calls from dPdTtr(). Thanks to Tomas Kalibera for the
      suggestion.

  }
}

\section{Changes in CHNOSZ version 1.0.7 (2015-11-19)}{
  \itemize{

    \item Add functions usrfig() (get figure limits in user coordinates) and
      label.figure() (add label to figure outside of plot region).

    \item Add demo copper.R (complexation of copper with glycine, uses
      mosaic()).

    \item Using new supporting function ibasis(), swap.basis() and mosaic()
      now allow basis species to be specified by name, not only formula.

    \item diagram() now has '...' argument to pass additional options to
      plot() (useful with diagram(tplot=FALSE, ...)).

    \item Update vignette 'equilibrium.Rnw' with better definitions of concepts,
      organization of functions, and examples and applications. Now uses
      knitr.

    \item Remove demo diagram.R (concepts better shown in equilibrium.Rnw).

    \item Rename basis.matrix() to basis.elements().

    \item Argument of species.basis() defaults to current species indices.
      ( species.basis() %*% basis.elements() returns elemental composition
      of currently defined species.)

    \item Move demos phosphate.R and nonideal.R to help page nonideal.Rd.

    \item Add options to thermo$opt: 'ideal.H' and 'ideal.e'. Default TRUE tells
      nonideal() to set activity coefficients of proton and electron to
      zero.

    \item Add 'method' argument to equilibrate() (can be 'reaction' or
      'boltzmann').

    \item For compatibility with a change in R-devel ("aperm(a, *) now preserves
      names(dim(a))"), modify slice.affinity() to preserve names of
      dimensions, and use check.attributes=FALSE or check.names=FALSE in
      some tests.

  }
}

\section{Changes in CHNOSZ version 1.0.6 (2015-10-19)}{
  \itemize{

    \item Rewrite rho.IAPWS95() to be able to invert density from IAPWS-95
      equations for a more extensive range of T,P values.

    \item Update R dependency to R-3.1.0, needed for 'extendInt' argument
      of uniroot() (used by rho.IAPWS95()).

    \item Add demo/density.R to show density of H2O calculated using
      rho.IAPWS95() (optionally using water.SUPCRT92()).

    \item Remove the relatively trivial demo/cordierite.R.

    \item Update wrapper function for parallel calculations palply() (export of
      variables, used in demo/density.R) and add tests in
      test-util.program.R.

    \item Tests that initiate calls to palply() (and therefore makeCluster())
      failed with 'cannot open file 'startup.Rs': No such file or directory'
      Fixed by adding Sys.setenv("R_TESTS" = "") to test-all.R.
      Issue discussed here:
      https://github.com/hadley/testthat/issues/129
      https://github.com/hadley/testthat/issues/144

    \item Fix bugs in which.pmax() that prevented proper assignment of
      attributes in output, and functionality for lists of length 1.

    \item mosaic() now multiplies affinities by _relative_ abundances of basis
      species when blend=TRUE. Also, fix bug where using blend=TRUE gave
      wrong results due to incorrect calculation of total activity of basis
      species. Thanks to Grayson Boyer for the bug reports.

    \item Add grDevices, graphics, stats, utils to Imports: in DESCRIPTION, and
      add specific importsFrom lines in NAMESPACE.

    \item Fix incorrect entry for entropy of aqueous methionine and [Met] in
      OBIGT.csv. Thanks to Apar Prasad for reporting this.

    \item Some fixes for compatibility with new version of testthat (increase
      tolerance for one test in test-eos.R because of corrected alignment
      actual/expected values in equals(); move inst/tests to
      inst/test/testthat; use test_check instead of test_package in
      test-all.R (for R CMD check)).

  }
}

\section{Changes in CHNOSZ version 1.0.5 (2015-05-19)}{
  \itemize{

    \item Following "Writing R Extensions", remove uncommented \code{\\donttest{}} code
      from Rd files. Some examples have had \code{\\donttest{}} markup removed
      without other modification, some non-representative examples have
      been deleted, and some longer examples have been moved to demo files.

    \item Some previous \code{\\donttest{}} examples have been moved to demos 'findit.R'
      (replaces old demo), and to new demos 'revisit.R', 'ionize.R',
      'buffer.R', 'wjd.R', 'yeastgfp.R', 'diagram.R'.

    \item Rd files now use standard R \code{\\keyword{}}s ('package' and 'utilities');
      non-standard keywords have been changed to \code{\\concept{}}s (see CHECKLIST
      for description, and CHNOSZ-package.Rd for overview).

  }
}

\section{Changes in CHNOSZ version 1.0.4 (2015-05-18)}{

  \subsection{NEW FEATURES}{
    \itemize{
      \item Add mosaic() function for affinity calculations with changing basis
        species. See also new 'mosaic.R' demo.

      \item Add thermo$opt$varP option for subcrt() to calculate Gibbs energies
        of gases using a variable-pressure standard state. Thanks to Kirt
        Robinson for the submitting the problem that led to this change.

      \item Add 'solubility.R' demo, and add calcite solubility plot to
        anintro.Rnw.
    }
  }

  \subsection{DATA CONSISTENCY}{
    \itemize{
      \item Move aqueous Ti(OH)4 and H4SiO4 to OBIGT-2.csv.

      \item Add warning about data consistency and use of examples to
        CHNOSZ-package.Rd and anintro.Rnw.
    }
  }

  \subsection{OTHER CHANGES}{
    \itemize{
      \item When making filled color diagrams, diagram() calls image() with
        useRaster=TRUE. This speeds up plotting, eliminates artifacts (i.e.
        blank lines) at some resolutions, and greatly reduces the size
        of PDF files.

      \item diagram() returns 'lx', 'ly' and 'is' for calculated positions of
        labels and indices of those species in predominance diagrams,
        including when the labels aren't plotted (by setting names=FALSE).

      \item In energy(), add missing 'IS' argument in call to subcrt(). Thanks to
        Grayson Boyer for the bug report.

      \item In nonideal(), skip proton and electron. Add test for calculations
        at constant and variable IS to test-util.affinity.R.

      \item Modify hotspring.Rnw to use knitr, and add example calculations of
        relative abundances of microbial phyla from Dick and Shock, 2013.

      \item Add files with average amino acid compositions of proteins from Bison
        Pool hot spring grouped according to annotation keyword (DS11.csv)
        (moved here from protein.csv) or phylum assignment (DS13.csv).

      \item water.lines() gets 'O2state' argument to specify state of O2.

      \item Add 'as.residue' argument to diagram(), to divide reactions by
        balance coefficients without rescaling to whole species formulas
        (cf. 'normalize' argument which uses rescaling).

      \item In equilibrate(), change name of 'stay.normal' argument to
        'as.residue'.

      \item In test-diagram.R, test that same diagrams are produced using
        'normalize' in equilibrate() or diagram().

      \item In test-diagram.R, test that same diagrams are produced using
        'as.residue' in equilibrate() or diagram().

      \item Add 'tplot' argument to diagram() (change to FALSE to use R's default
        plot set-up).

      \item R version dependency changed to (>= 3.0.0).
    }
  }

}

\section{Changes in CHNOSZ version 1.0.3 (2014-01-12)}{
  \itemize{

    \item Updated extdata/protein/Sce.csv.xz using Saccharomyces Genome
      Database protein_properties.tab and SGD_features.tab dated
      2013-08-24.

    \item more.aa() includes SGDID and gene name in the abbrv and ref columns.

    \item In extdata/refseq, scripts and data files were updated for NCBI
      Reference Sequence (RefSeq) release 61 (2013-09-09). Code was
      adapted to deal with WP multispecies accessions.

    \item read.fasta() gets new argument 'id'; when supplied, it is used for
      the protein names in the output, in place of those read from the
      FASTA headers.

    \item When reading protein names from the FASTA headers, read.fasta()
      stops only at the first space, not space or underscore as before.

    \item info() no longer specially sets state of "O2" to gas. The name
      "oxygen", or the combination ("O2", "gas"), can be used to retrieve
      data for the gas.

    \item In thermo$obigt, names of gases (e.g. "oxygen") are used only for
      the gaseous species; names were removed from dissolved species in
      rows 63-79. Both gases and aqueous species can continue to be
      referenced by their chemical formulas.

    \item In read.expr(), allow multiple filter specifications.

    \item revisit() has new argument loga0, a single vector of base-10
      logarithms of activities of species used to calculate the base-2 log
      ratio ( log2(a1/a0) ).

    \item Updated tests to be compatible with testthat version 0.8 (in
      development). Some tests now check for e.g. an error *and* warning
      produced by a single function call.

  }
}

\section{Changes in CHNOSZ version 1.0.1 (2013-07-04)}{
  \itemize{

    \item Allow IAPWS-95 calculations to be used in subcrt(): in 
      water.IAPWS95(), rename internal function from 'epsilon' to 'diel';
      in water(), return upper-case versions of names of properties.

    \item rho.IAPWS95() has modified search interval to give results at
      T < 250 K and P < 300 bar. Thanks to Marc Neveu.

    \item In water.AW90(), analytical solution to polynomial equation is
      now used to calculate dielectric constant. Thanks to Marc Neveu.

    \item count.aa() and read.fasta() handle DNA sequences, with new argument
      'type' (protein or DNA).

    \item nucleicacids() split into nucleic.formula() and nucleic.complement().

  }
}

\section{Changes in CHNOSZ version 1.0.0 (2013-03-28)}{

  \subsection{MAJOR USER-VISIBLE CHANGE}{
    \itemize{
      \item On attaching the package, the user is asked to load the 'thermo'
        object, containing thermodynamic data and system settings, using
        data(thermo). This is required because packages are not permitted to
        alter the search path (but the user may).
    }
  }

  \subsection{OTHER CHANGES}{
    \itemize{
      \item Fix calculation of free energy derivative in wjd().

      \item Add 'stay.normal' argument to equilibrate().

      \item If obigt$G is available, hkf() returns this value, not NA, at Tr, Pr.

      \item mod.obigt() defaults to taking chemical formula from the species name,
        and checks for validity of formula.

      \item Add example for LYSC_CHICK to protein.info.Rd.

      \item Enable DGtr in findit(), 1-D DGtr plot in revisit().

      \item Disable a check in valTP function of H2O92D.f to allow properties of
        H2O to be calculated below 0.01 degrees C, to -20 degrees C.

      \item Remove thermo$water and thermo$water2 storage of previous results;
        they gave no significant speed gain in running examples.

      \item Split IAPWS95.Rd out of water.Rd.

      \item expr.* functions now do not wrap their values in as.expression()
        (makes compounding expressions, e.g. with substitute(), easier).

      \item Rename dl.aa() to uniprot.aa().

      \item Add start and stop arguments to count.aa(), read.fasta(),
        uniprot.aa().

      \item add.protein() replaces amino acid compositions for existing
        proteins with the same name.

      \item Examples of calculation of affinity of formation of CSG_METVO
        (following Dick and Shock, 2011) added to protein.info.Rd.

      \item Use consistent names for water properties (Speed, diel, QBorn, ...).

      \item Add water.props() to get names of properties of water.

      \item Remove 'isat' argument from water.SUPCRT92(); function now accepts
        'Psat' as value for 'P' argument.

      \item Separate rho.IAPWS95() from water.IAPWS95().

      \item In addition to the original regression variables, EOSvar() recognizes
        names for available properties in water(), or can use the name to
        get a user-defined function of temperature and pressure.

      \item Simplify EOSplot() somewhat (don't group data by pressure ranges).

      \item EOSlab() gets label from attribute (if present) of user-defined
        function.

      \item Modify test-util.data.R and wjd.Rd to pass R CMD check using R
        configured with --disable-long-double .

      \item The name of the environment affected by data(thermo), and used in many
        functions, is changed from CHNOSZ:thermo to CHNOSZ.

      \item Remove read.supcrt() and write.supcrt().

      \item guess() now defaults to "stoich" method, not "central". (Needed since
        limSolve package is not found during R-Forge checks on Windows.)
    }
  }

}

\section{Changes in CHNOSZ version 0.9-9 (2013-01-01)}{

  \subsection{MAJOR CHANGES}{
    \itemize{
      \item Split the functionality of diagram() into equilibrate() and diagram().
        Old workflow: a <- affinity(); d <- diagram(a) .
        New workflow: a <- affinity(); e <- equilibrate(a); d <- diagram(e) .
        Old workflow is still usable for plotting the values of affinity, or
        for making predominance diagrams using the maximum affinity method.

      \item The 'thermo' object, which holds the thermodynamic database,
        and system definitions (made by the user), is now placed in an
        environment named 'CHNOSZ:thermo' on loading the package. Therefore, 
        instances of '<<-' in the code now refer to this environment
        instead of the global environment.

      \item Create a set of tests in 'inst/tests', particularly for functions that
        have been modified during this development cycle, and add a Suggests 
        dependency on 'testthat'.

      \item Move the code for the temperature and pressure derivatives of the "g"
        function (related to the solvation parameter omega) to a new function 
        gfun(); incorporate some fixes and a series of test_that() tests.
        There is some impact on the calculated Gibbs energies of charged
        species.
    }
  }

  \subsection{NEW FUNCTIONS}{
    \itemize{
      \item Add wjd() implementing the steepest descent algorithm for free energy
        minimization described by White et al., 1958. Also add supporting 
        functions element.potentials(), invertible.combs() for finding 
        linearly independent combinations of rows of a matrix, 
        is.near.equil(), and run.wjd().

      \item Add guess() as another supporting function for wjd(), to produce 
        initial guesses of moles of species satisfying a given elemental
        bulk composition, and a Suggests dependency on 'limSolve'.

      \item New function i2A() for generating a stoichiometric matrix from indices
        of species in the thermodynamic database.

      \item Add protein.equil() for step-by-step calculation of chemical
        activities of proteins in metastable equilibrium.

      \item Add an objective function DGtr() for calculating the Gibbs energy of
        transformation of a system at constant temperature, pressure and
        chemical activities of basis species.

      \item Add msgout(), which is a modification of message() from base R. Now
        used instead of cat(), to allow suppressing messages (e.g. during
        testing with test_that).
    }
  }

  \subsection{CHANGES TO ARGUMENTS}{
    \itemize{
      \item In equilibrate(), the argument 'logact' (inherited from diagram()),
        specifying the logarithm of activity of the balanced quantity, has
        been renamed to 'loga.balance'. In the result, rename 'logact' to
        'loga.equil', containing the equilibrium logarithms of activities
        of the species of interest.

      \item The 'residue' argument of diagram() has been changed to 'normalize'
        in diagram() and equilibrate(). normalize=FALSE is always the default,
        including for systems of proteins.

      \item In balance(), the value of 'balance' used to indicate protein length
        has been changed from 'PBB' to 'length'.

      \item Everywhere it appeared, the logical argument 'do.plot' has been 
        renamed to 'plot.it' (diagram(), revisit(), findit(), transfer()).
        (This scheme is more consistent with e.g. qqplot().)

      \item 'do.phases' argument in subcrt() and affinity() has been renamed to
        'exceed.Ttr'. When that argument is FALSE (the default for subcrt()),
        the Gibbs energies of mineral phases at temperatures beyond their 
        transition temperature are set to NA, instead of 999999 used 
        previously.

      \item ZC() now accepts a numeric argument, referring to one or more species
        indices in the thermodynamic database.
    }
  }

  \subsection{CHANGES TO OUTPUT}{
    \itemize{
      \item subcrt() now only outputs T (temperature), P (pressure) and rho
        (density) columns if there is more than one T-P point ... makes
        unlist()ing the results easier (used in element.mu() and 
        basis.logact()).

      \item subcrt() now outputs NA values for properties at temperatures above
        the critical temperature of H2O, when Psat is being used.

      \item read.blast() accepts NA for 'similarity', 'evalue' and 'max.hits' 
        options. Descriptive column names are now assigned to the data frame
        returned by the function.

      \item energy.args() (called by affinity()) shows units in messages about
        limits of variables.

      \item EOSvar() has new variables invPPsi and invPPsiTTheta, used for
        temperature- and pressure-dependent regressions in the revised HKF
        equations of state.

      \item thermo.plot.new() now saves the graphics device parameters 
        (par(no.readonly=TRUE)) to thermo$opar the first time the function is 
        called, allowing the parameters to be restored after running examples
        that change them.

      \item Rewrite mod.obigt() (it's now used by info() when adding proteins) and
        add today() for returning today's date in the format used in SUPCRT 
        files.
    }
  }

  \subsection{REFACTORING OF FUNCTIONS}{
    \itemize{
      \item Split the primary functionality of makeup(), parsing of chemical 
        formulas, into a smaller makeup() function and supporting functions 
        count.charge(), count.formulas(), and count.elements(). The new 
        functions make extensive use of regular expressions, and no data 
        frames. Running makeup() over the ca. 3000 formulas in thermo$obigt 
        drops from ~35 to ~5 seconds on one machine.

      \item New function as.chemical.formula() to replace the previous
        functionality of makeup() for making string representations of
        chemical formulas.

      \item Replace element() with two separate functions, mass() and entropy().

      \item Units setting interface is now split between three separate functions:
        P.units(), T.units() and E.units().

      \item Replace describe() with describe.basis(), describe.property(),
        describe.reaction(). It is now fairly easy to make legends showing
        temperature, pressure and chemical activities with italic symbols,
        subscripts, and units.

      \item Reorganize axis.label() and its supporting functions 
        (now expr.species(), expr.property() and expr.units()). Add an example
        showing a plot annotated with chemical formulas and reactions.

      \item The functionality of info() is split into info(), info.character(),
        info.approx(), info.numeric() and info.text(). For ease of use, single
        approximate matches are accepted by info(), and searches for 'H2O' in
        the 'aq' state now return H2O(liq).

      \item The basis definition functionality of basis() is split into basis(), 
        is.basis(), put.basis(), mod.basis(), preset.basis() and 
        preset.logact().

      \item The basis swapping functionality of basis() is split into 
        basis.matrix(), element.mu(), basis.logact() and swap.basis(). Rename
        basis.comp() to species.basis(), and remove expand.formulas().

      \item Split aminoacids() into aminoacids() and count.aa().

      \item The monolithic protein() function no longer exists; it has been
        superseded by iprotein(), ip2aa(), get.aa(), dl.aa(), read.aa(),
        sum.aa(), and aa2eos(). The user shouldn't notice significant changes
        (other than in the composition of messages) when including proteins in
        functions like subcrt() and affinity().

      \item Functionality of get.protein() is split into more.aa() (amino acid
        compositions from model organisms) and stress() (proteins identified
        in stress response experiments).

      \item get.expr(), for reading abundances or expression levels of proteins
        from variously formatted data files, is renamed to read.expr(), with
        retrieval of amino acid compositions of proteins moved to more.aa().

      \item Replace ionize() with a completely rewritten and much easier to use
        ionize.aa() for calculating the additive ionization properties of 
        proteins; also add A.ionization(), usually invoked by affinity()

      \item residue.info() mostly replaced by protein.basis().

      \item Rename equil.react() to equil.reaction(), and give it parallel
        potential (via palply) and a more efficient algorithm for determining
        limits of the uniroot search. Also remove a redundant argument; 
        arguments are now identical to those of equil.boltzmann() (and that
        function was renamed from equil.boltz()).

      \item Objective functions used in revisit() and findit() now each have
        their own definitions, with an attribute indicating whether the
        function is minimized or maximized.

      \item Separate idealgas.IAPWS95() and residual.IAPWS95() from
        water.IAPWS95() (in order to write some test_that tests.)
    }
  }

  \subsection{DATA UPDATES}{
    \itemize{
      \item Add 148 liquid and 148 crystalline acyclic isoprenoids, polycyclic 
        alkanes, polynuclear aromatic hydrocarbons (PAH) and 62 crystalline
        double ether-bonded or ester-bonded n-alkanes from Tables 16-23 of 
        Richard and Helgeson, 1998. The properties of the following 
        crystalline and liquid compounds, taken from those tables, replace the
        previously entered values from the "reference model compound" tables
        appearing earlier in that paper: biphenyl, naphthalene, 
        1-methylnaphthalene, 1,8-dimethylnaphthalene, 2,3-dimethylnaphthalene
        (cr only), anthracene, pyrene.

      \item Add 13 crystalline, 29 liquid and 39 gaseous organic iodine compounds
        from Richard and Gaona, 2011.

      \item Add crystalline peptide sidechain and [AABB] and [PBB] groups and 
        dipeptide model compounds and revised equations of state parameters 
        for crystalline leucine and revised standard Gibbs energies of
        crystalline and aqueous methionine and its aqueous sidechain group
        [Met], from LaRowe and Dick, 2012. Also, use revised volumes of all
        other crystalline amino acids given in that paper together with those
        of Lys:HCl and His:HCl calculated using an effective volume of HCl(cr)
        equal to 33.8 cm3/mol (difference between Arg and Arg:HCl).

      \item Add 6 aqueous chloroethylenes from Haas and Shock, 1999.

      \item Move superseded aqueous methionine sidechain group [Met] to 
        OBIGT-2.csv so that its properties are available for reproduction
        of published results (relevant to some examples in the package).

      \item As revised makeup() function now strictly interprets a signed value at
        the end of a chemical formula as charge, the formula of the electron
        in thermo$obigt is changed from "Z-1" to "Z0-1".

      \item In thermo$obigt, change the value for standard Gibbs energy of H2O 
        from -56688 to -56687.711481 cal mol-1, to be consistent with the 
        value generated by the fortran code from SUPCRT92 (H2O92D.f). Although
        the latter is the source of the properties of H2O for many functions
        in the package, there is an occasional function that accesses the
        tabulated value at 25 degrees C (e.g., element.mu() and 
        basis.logact()).

      \item In protein.csv, change organism code BACST (Bacillus 
        stearothermophilus) to GEOSE (Geobacillus stearothermophilus) for SLAP
        and DPO1 proteins, and also apply changes in vignettes and examples.

      \item In stress.csv, change ECO to Eco and SGD to Sce.
    }
  }

  \subsection{EXTDATA UPDATES}{
    \itemize{
      \item In extdata/refseq, scripts and data files were updated for NCBI
        Reference Sequence (RefSeq) release 55 (2012-09-17).

      \item In extdata/bison, sample BLAST output files for Bison Pool metagenome
        use target database generated from RefSeq release 55.

      \item Add P(ressure) column to extdata/cpetc/SOJSH.csv and a stopifnot() 
        test for similarity to the experimental data to the example in 
        water.Rd.

      \item In extdata/protein rename ECO.csv.xz to Eco.csv.xz and SGD.csv.xz to
        Sce.csv.xz.

      \item In extdata/thermo add RH98_Table15.csv; this file together with new
        function RH2obigt() is used to calculate thermodynamic properties of
        organic compounds using group contributions from Richard and 
        Helgeson, 1998.
    }
  }

  \subsection{DOCUMENTATION AND VIGNETTES}{
    \itemize{
      \item Add a Known Bugs section to the package help page.

      \item Move 'extra' examples previously available in longex() to individual
        demos, and add a function demos() to run them all.

      \item All examples now run without any warnings (at least, as intended).

      \item Add an example based on Shock and Canovas, 2010 for the 'transect' 
        mode of affinity().

      \item ionize.aa() has examples of contour plots as a function of temperature
        and pH.

      \item Group residue.formula(), protein.formula(), protein.name(), 
        protein.length() into new documentation topic (util.protein.Rd).

      \item Add 'sideeffects.Rd' to document some of the side effects of 
        functions.

      \item Add 'objective.Rd' for the objective functions.

      \item To make output of examples reproducible, change mod.obigt() and 
        protein() to not include current date/time for new data entries 
        (remains an option for mod.obigt()).

      \item Move vignette sources to 'vignettes' directory.

      \item Remove vignettes formation.Rnw and xadditivity.Rnw 
        (and, for the latter, the package Suggests dependency on xtable).

      \item Rename protactive.Rnw to equilibrium.Rnw, with some changes.

      \item Add vignette wjd.Rnw to accompany the new function wjd().

      \item Add a release CHECKLIST to the installation directory of the package.

      \item Replace "degrees" with UTF-8 degree symbol in second argument of
        \code{\\eqn{}{}}; add \code{\\encoding{UTF-8}} to affected Rd files.
    }
  }

  \subsection{OTHER CHANGES}{
    \itemize{
      \item A matrix is now returned by makeup(), GHS() and some other functions
        to avoid the performace penalty associated with data frames.

      \item Add parallel as a Suggests dependency and replace mylapply() (based on
        'multicore') with palply(), as a wrapper for parLapply and lapply.
        palply() now invokes parLapply for length(X) > 10000 instead of 100.

      \item Comment out lines containing WRITE and STOP statements in src/H2O92D.f
        since they are discouraged by CRAN guidelines (and calls to them are
        unlikely to be encountered while using the package).

      \item Remove thermo$opt$verbose and thermo$opt$online.

      \item In c2s(), use 'collapse' argument of paste() instead of a for loop.

      \item Change R dependency to R >= 2.12.0, required for useDynLib in
        NAMESPACE to find the shared library on Windows.
    }
  }

}

\section{Changes in CHNOSZ version 0.9-7 (2011-08-23)}{
  \itemize{

    \item Restore ZipData: no in DESCRIPTION, needed for building the package
      on R < 2.13.0 on Windows.

    \item Remove some incorrect titles for barplots in diagram().

    \item Add more informative titles for plots in first example of diagram.Rd.

    \item In EOSvar() and EOSlab(), name the first argument 'EXPR' in calls to 
      switch() to avoid partial matching with 'E' (expansivity).

    \item Shorten CHNOSZ-package.Rd and clean up formatting of some examples
      in other documentation topics.

  }
}

\section{Changes in CHNOSZ version 0.9-6 (2011-08-18)}{

  \subsection{SIGNIFICANT CHANGES TO FUNCTIONS}{
    \itemize{
      \item New function browse.refs() for listing references for thermodynamic
        data and opening associated URLs. Supported by extdata/js/sorttable.js 
        for sorting table of references in browser.

      \item Add EOSregress(), EOSvar(), EOSlab(), EOSplot(), EOScoeffs() for 
        regressing and comparing equations-of-state parameters from 
        heat capacity and volume data for aqueous species.

      \item Add anim.plasma() and anim.carboxylase() for making animated GIFs of
        equilibrium activity diagrams. All anim.* functions now also work
        on Windows (with ImageMagick installed).

      \item New functions checkEOS(), checkGHS() for checking self-consistency of
        individual database entries. check.obigt() replaces 'check' argument 
        of info().

      \item New argument 'chains' for protein() can be used to specify the number
        of polypeptide chains in group additivity calculations.

      \item protein() now looks for protein backbone group named '[PBB]' for group
        additivity of crystalline proteins.

      \item Remove count.taxa(), splitting its functionality into read.blast()
        and id.blast(); new function write.blast().

      \item eqdata() can now also extract results for solids, saturation states of
        minerals, and major speciation of basis species from EQ6 output files.
    }
  }

  \subsection{DATA UPDATES}{
    \itemize{
      \item Rename 'source.csv' to 'refs.csv'. Add URLs for most references.

      \item Column names 'source1' and 'source2' in OBIGT.csv, changed to
        'ref1' and 'ref2'.

      \item In protein.csv, add amino acid compositions of 105 
        model proteins derived from metagenomic sequences at 
        Bison Pool (Dick and Shock, 2011).
        protein: overall, transferase, transport, dehydrogenase, ...;
        organism: bisonN, bisonS, bisonR, bisonQ, bisonP

      \item Add properties for aqueous bromine (Br2) and iodine (I2) from 
        Wagman et al., 1982.

      \item Add properties of crystalline and liquid groups from
        Helgeson et al., 1998 and Richard and Helgeson, 1998.
    }
  }

  \subsection{EXTDATA UPDATES}{
    \itemize{
      \item Reorganize the 'extdata' directory, putting all files into one
        of nine subdirectories (abundance, bison, cpetc, fasta, js, protein,
        refseq, taxonomy, thermo).

      \item To reduce installed package size, compress many data files in extdata 
        using xz.

      \item In extdata/bison, update partial protein BLAST output files for five 
        Bison Pool sites to use target database made up of microbial protein 
        sequences in NCBI Reference Sequence database version 47.

      \item In extdata/cpetc, new files 'Cp.CH4.HW97.csv' and 'V.CH4.HWM96.csv'
        contain experimental heat capacity and volume data for aqueous methane
        from Hnedkovsky and Wood, 1997 and Hnedkovsky et al., 1996.

      \item In extdata/refseq, rename 'taxid_phylum.csv' to 'taxid_names.csv'.

      \item In extdata/thermo/OBIGT-2.csv, fix incorrectly entered values of c1 from 
        Dalla-Betta and Schulte, 2010. Also apply corrected values of c1 and c2 
        and minor corrections for delta-H and S from Marini and Accornero, 2010.
    }
  }

  \subsection{DOCUMENTATION AND VIGNETTES}{
    \itemize{
      \item Change vignette file names and titles to be in the same alphabetical
        order to achieve the desired sorting in browseVignettes().

      \item Add .Rinstignore to inst/doc to exclude figure, bibliography and
        LyX files from the installation directory.

      \item Adapted 'formation.Rnw' to work with pdflatex, and remove
        Makefile in inst/doc.

      \item Add 'hs-chemistry.Rnw' with calculations of relative stabilities 
        of model proteins in a hot spring (Dick and Shock, 2011).

      \item Add 'protactiv.Rnw' with detailed information on calculating equilibrium
        activities of proteins, and comparisons with experimental abundances in 
        human blood plasma and E. coli.

      \item Add 'xadditivity.Rnw' with examples of group contribution calculations
        (using qr.solve()) and comparison of two group contribution schemes.
    }
  }

  \subsection{OTHER CHANGES}{
    \itemize{
      \item In DESCRIPTION, change Depends to a more recent R (>= 2.10.0),
        required for reading of compressed (including .xz) files.

      \item Add Suggests: xtable in DESCRIPTION (used for 'xadditivity' vignette).

      \item In DESCRIPTION, change "ZipData" to "BuildResaveData", keeping 
        the setting at "no". Resaving data as compressed files reduces the
        installed package size by about 200Kb, but complicates showing the file 
        formats to new users.

      \item Loading of shared object for calculating properties of H2O
        now done by useDynLib in NAMESPACE.

      \item Change .First.lib() to .onAttach(). In .onAttach, use 
        packageStartupMessage() instead of cat().
    }
  }

}
