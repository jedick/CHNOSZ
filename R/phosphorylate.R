# Calculate affinity of phosphorylation reactions taking account of speciation
# 20251206 first version (extracted from sugars paper script) jmd
# 20251208 add phospho_plot()

phosphorylate <- function(reactant, P_source, loga_reactant = 0, loga_product = 0, loga_P_source = 0, loga_P_remainder = 0, const_pH = 7, ...) {

  # Affinity is calculated by summing:
  # 1) reactant + P = product + H2O (m1) 
  # 2) P_source + H2O = P_remainder + P (m2)
  # NOTE: reaction 2 is defined in reverse (to form H2O), then the opposite of the affinity is used for the sum

  # P_source can be "P" (basic reaction) or "PP", "acetylphosphate", or "ATP" (extended reaction)
  # Mapping from P_source -> P_remainder:
  # P (H3PO4) -> H2O
  # PP (H4P2O7) -> H3PO4
  # acetylphosphate -> acetic acid
  # ATP -> AMP
  # NOTE: loga_P_remainder applies to all of these *except* H2O

  # Terminology:
  # Complex species are species that take multiple forms (e.g. ionization or oxidation states) depending on pH or other variables
  # A basic reaction has complex species that are valid basis species (independent components)
  # --> This is what mosaic() deals with
  # An extended reaction has complex species all of which cannot be used as basis species (they are non-independent)
  # However, the energetics of an extended reaction can be modeled as a sum of basic reactions
  # --> This is what this function deals with

  # Examples:
  # Complex species:    acetic acid = [acetic acid + acetate]
  # Basic reaction:     acetic acid + P = acetylphosphate + H2O
  # Extended reaction:  acetic acid + PP = acetylphosphate + P

  ## Setup reaction 1
  if(reactant == "acetic acid") {
    # Basic reaction: acetic acid + P = acetylphosphate + H2O
    # Load initial species for mosaic reaction (uncharged species)
    basis(c("acetic acid", "H3PO4", "acetylphosphate0", "O2", "H+"))
    # The basis species we will speciate using mosaic()
    bases <- list(
      c("H3PO4", "H2PO4-", "HPO4-2", "PO4-3"),
      c("acetic acid", "acetate"),
      c("acetylphosphate0", "acetylphosphate-1", "acetylphosphate-2", "acetylphosphate-3")
    )
  } else if(reactant == "glycerol") {
    # Basic reaction: glycerol + P = 1-glycerolphosphate + H2O
    basis(c("glycerol", "H3PO4", "3-phosphoglycerol", "O2", "H+"))
    bases <- list(
      c("H3PO4", "H2PO4-", "HPO4-2", "PO4-3"),
      c("3-phosphoglycerol", "3-phosphoglycerol-1", "3-phosphoglycerol-2")
    )
  } else if(reactant == "adenosine_to_AMP") {
    # Basic reaction: adenosine + P = AMP + H2O
    basis(c("adenosine", "H3PO4", "H2AMP", "N2", "O2", "H+"))
    bases <- list(
      c("H3PO4", "H2PO4-", "HPO4-2", "PO4-3"),
      c("H2AMP", "HAMP-", "AMP-2")
    )
  } else if(reactant == "adenosine_for_RNA") {
    # Basic reaction: adenosine + P = AMP + H2O
    # To reproduce calculations from LaRowe and Dick (2025):
    # PO4-3 is not included because monophosphate nucleotides can't have a -3 charge 20250426
    # AMP-2 is not included because the repeating unit in RNA can't have this charge 20250424
    basis(c("adenosine", "H3PO4", "H2AMP", "N2", "O2", "H+"))
    bases <- list(
      c("H3PO4", "H2PO4-", "HPO4-2"),
      c("H2AMP", "HAMP-")
    )
  } else if(reactant == "adenosine_to_cAMP") {
    # Basic reaction: adenosine + P = cAMP + H2O
    basis(c("adenosine", "H3PO4", "cyclic-HAMP", "N2", "O2", "H+"))
    bases <- list(
      c("H3PO4", "H2PO4-", "HPO4-2", "PO4-3"),
      c("cyclic-HAMP", "cyclic-AMP-1")
    )
  } else if(reactant == "ribose") {
    # Basic reaction: ribose + P = ribose-5-phosphate + H2O
    basis(c("ribose", "H3PO4", "ribose-5-phosphate", "O2", "H+"))
    bases <- list(
      c("H3PO4", "H2PO4-", "HPO4-2", "PO4-3"),
      c("ribose-5-phosphate", "ribose-5-phosphate-1", "ribose-5-phosphate-2")
    )
  } else if(reactant == "uridine") {
    # Basic reaction: uridine + P = UMP + H2O
    basis(c("uridine", "H3PO4", "H2UMP", "N2", "O2", "H+"))
    bases <- list(
      c("H3PO4", "H2PO4-", "HPO4-2", "PO4-3"),
      c("H2UMP", "HUMP-", "UMP-2")
    )
  } else if(reactant == "AMP") {
    # Basic reaction: AMP + P = ADP + H2O
    basis(c("H2AMP", "H3PO4", "H3ADP", "N2", "O2", "H+"))
    bases <- list(
      c("H3PO4", "H2PO4-", "HPO4-2", "PO4-3"),
      c("H2AMP", "HAMP-", "AMP-2"),
      c("H3ADP", "H2ADP-", "HADP-2", "ADP-3")
    )
  } else if(reactant == "ADP") {
    # Basic reaction: ADP + P = ATP + H2O
    basis(c("H3ADP", "H3PO4", "H4ATP", "N2", "O2", "H+"))
    bases <- list(
      c("H3PO4", "H2PO4-", "HPO4-2", "PO4-3"),
      c("H3ADP", "H2ADP-", "HADP-2", "ADP-3"),
      c("H4ATP", "H3ATP-", "H2ATP-2", "HATP-3", "ATP-4")
    )
  } else if(reactant == "glucose") {
    # Basic reaction: glucose + P = glucose-6-phosphate + H2O
    basis(c("glucose", "H3PO4", "glucose-6-phosphate", "O2", "H+"))
    bases <- list(
      c("H3PO4", "H2PO4-", "HPO4-2", "PO4-3"),
      c("glucose-6-phosphate", "glucose-6-phosphate-1", "glucose-6-phosphate-2")
    )
  } else if(reactant == "pyruvic acid") {
    # Basic reaction: pyruvic acid + P = phosphoenolpyruvate + H2O
    basis(c("pyruvic acid", "H3PO4", "phosphoenolpyruvate", "O2", "H+"))
    bases <- list(
      c("H3PO4", "H2PO4-", "HPO4-2", "PO4-3"),
      c("pyruvic acid", "pyruvate"),
      c("phosphoenolpyruvate", "phosphoenolpyruvate-1", "phosphoenolpyruvate-2", "phosphoenolpyruvate-3")
    )
  } else {
    stop(paste("unrecognized reactant:", reactant))
  }

  # Make sure we are speciating H3PO4 (in addition to the other reactants)
  stopifnot(bases[[1]][1] == "H3PO4")

  # Set activity of H3PO4
  # Basic reaction: H3PO4 is the P_source
  if(P_source == "P") basis("H3PO4", loga_P_source)
  # Extended reaction with PP: H3PO4 is the P_remainder
  if(P_source == "PP") basis("H3PO4", loga_P_remainder)
  # Other extended reactions [don't set the H3PO4 activity]
  #   There is no H3PO4 in the overall reaction, but is present in the half reactions.
  #   For correct cancellation, we just use the default activity of H3PO4.

  # NOTE: The initial basis definition uses the name, but setting the activity uses the elemental formula
  formula_product <- rownames(basis())[3]
  basis(formula_product, loga_product)
  basis("pH", const_pH)
  # This is used for changing the activity of the reactant
  formula_reactant <- rownames(basis())[1]
  basis(formula_reactant, loga_reactant)

  # Generate overall reaction
  # Don't use e.g. species("acetylphosphate0"), because that reaction is just acetylphosphate0 = acetylphosphate0
  # Instead, form H2O from the basis species to generate the overall reaction
  species("H2O")

  # Calculate affinity of forming product from predominant basis species
  m1 <- mosaic(bases, ...)
  # For cAMP, double the reaction to consume one H3PO4 20251201
  if(reactant == "adenosine_to_cAMP") m1$A.species$values <- lapply(m1$A.species$values, "*", 2)

  ## Setup reaction 2
  if(P_source == "P") {
    # We just want the first reaction (the basic reaction)
    m2 <- NULL
    a12 <- m1$A.species$values[[1]]
    P_reaction <- NULL
  } else {
    # Save the existing basis and species definition (used for m1)
    basis_1 <- basis()
    species_1 <- species()
    # Mosaic for opposite of reaction 2:
    if(P_source == "PP") {
      # Form H2O from H4P2O7 and H3PO4
      basis(c("H4P2O7", "H3PO4", "O2", "H+"))
      basis("H4P2O7", loga_P_source)
      basis("H3PO4", loga_P_remainder)
      bases <- list(
        c("H4P2O7", "H3P2O7-", "H2P2O7-2", "HP2O7-3", "P2O7-4"),
        c("H3PO4", "H2PO4-", "HPO4-2", "PO4-3")
      )
    } else if(P_source == "acetylphosphate") {
      # Form H2O from acetylphosphate and acetic acid
      basis(c("acetylphosphate0", "acetic acid", "H3PO4", "O2", "H+"))
      basis("C2H5O5P", loga_P_source)
      basis("C2H4O2", loga_P_remainder)
      bases <- list(
        c("acetylphosphate0", "acetylphosphate-1", "acetylphosphate-2", "acetylphosphate-3"),
        c("acetic acid", "acetate"),
        c("H3PO4", "H2PO4-", "HPO4-2", "PO4-3")
      )
    } else if(P_source == "ATP") {
      # Form H2O from ADP and ATP
      basis(c("H4ATP", "H3ADP", "H3PO4", "N2", "O2", "H+"))
      basis("C10H16N5O13P3", loga_P_source)
      basis("C10H15N5O10P2", loga_P_remainder)
      bases <- list(
        c("H4ATP", "H3ATP-", "H2ATP-2", "HATP-3", "ATP-4"),
        c("H3ADP", "H2ADP-", "HADP-2", "ADP-3"),
        c("H3PO4", "H2PO4-", "HPO4-2", "PO4-3")
      )
    }
    species("H2O")
    basis("pH", const_pH)
    # Use the bases defined above and take the other arguments (e.g. pH, T, and P) from m1
    m_args <- m1$args
    m_args[["bases"]] <- bases
    # TODO: handle ionic strength
    #if(IS > 0) m_args[["IS"]] <- IS
    m2 <- do.call(mosaic, m_args)
    # The affinity of the extended reaction
    a12 <- m1$A.species$values[[1]] - m2$A.species$values[[1]]
    # Return the formation reaction to use in describe.reaction() 20251201
    P_reaction <- species()
    # Restore the previous basis and species definitions for later calculations
    thermo(basis = basis_1)
    thermo(species = species_1)
  }

  # Put together the output
  list(m1 = m1, m2 = m2, a12 = a12, P_reaction = P_reaction)

} # end of phosphorylate()

# Define a function to make the plots for a given reaction
phospho_plot <- function(reactant, P_source) {

  # Reaction-independent settings (activities of species)
  # The product (phosphorylated species)
  loga_product <- -6
  # The reactant: we will loop over these to make a series of plots
  logas_reactant <- c(-6, -4, -2)

  # For the RNA model described in LaRowe and Dick (2025), loga_P_source (H3PO4) is set to loga_reactant
  # Otherwise, use constant values for loga_P_source and loga_P_remainder
  if(reactant != "adenosine_for_RNA") {
    loga_P_source <- -3
    loga_P_remainder <- -3
  }

  # Plot settings
  pH <- c(0, 12)
  T <- c(0, 200)
  P <- c(1, 5000)
  res <- 50

  # Define which contour levels to show
  levels <- c(-1e6, seq(-100, 100, 10), 1e6)
  # Use thick line for 0
  lwd <- ifelse(levels == 0, 2, 1)
  # Use shades of blue for exergonic, white for endergonic
  # Take away the 2 lightest and 1 darkest shades for better readability
  blues <- hcl.colors(14, "Blues")[2:12]
  col <- c(blues, rep("#FFFFFF", 11))

  # Use the top row (panel 7) for the reaction label
  top <- t(matrix(rep(7, 3)))
  bottom <- matrix(1:6, nrow = 2, byrow = TRUE)
  mat <- rbind(top, bottom)
  layout(mat, heights = c(1, 8, 8))

  # Loop over temperature and pressure for rows of figure
  for(yvar in c("T", "P")) {

    for(iloga in seq_along(logas_reactant)) {

      # Get single reactant loga
      loga_reactant <- logas_reactant[iloga]
      if(reactant == "adenosine_for_RNA") {
        # This is used to reproduce calculations from LaRowe and Dick (2025)
        loga_P_source <- loga_reactant
        # With H3PO4 as P_source there is no P_remainder, and this setting has no effect
        loga_P_remainder <- NA
      }

      # Perform calculations and define diagram settings for temperature or pressure
      if(yvar == "T") {
        # Calculate affinity of forming product from predominant basis species as a function of pH and temperature
        result <- phosphorylate(reactant, P_source, loga_reactant, loga_product, loga_P_source, loga_P_remainder, pH = c(pH, res), T = c(T, res))
        # Method for labeling contour lines
        method <- "flattest"
        # Legend placement, space, and expression
        legend.x <- "topright"
        legend.space <- "   "
        legend.expr <- as.expression(quote(italic(P)[SAT]))
        # Label offset (0 for a-c)
        dlab <- 0
      }
      if(yvar == "P") {
        # Use P_source=P_source to avoid argument collision with 'P' (pressure)
        result <- phosphorylate(reactant, P_source=P_source, loga_reactant, loga_product, loga_P_source, loga_P_remainder, pH = c(pH, res), P = c(P, res))
        method <- "edge"
        legend.x <- "bottomright"
        legend.space <- "      "
        legend.expr <- as.expression(quote(25~degree~C))
        # Label offset (3 for d-f)
        dlab <- 3
      }

      # Use basic mosaic output to get the plotting variables
      a <- result$m1$A.species
      # The next line is a workaround for y-axis mislabeled as log a ΣP
      # (sum of P activity in different basis species - but we want P(bar)) 20250422
      a$basis <- a$basis[colnames(a$basis) != "P"]
      # Start with blank diagram
      diagram(a, names = NA)

      # Get temperature values in Kelvin
      TK <- convert(a$vals$T, "K")
      if(yvar == "P") {
         # Label tick mark for 1 bar
         axis(2, at = 1)
         # Use constant T for yvar == "P"
         TK <- convert(25, "K")
      }

      # The affinity of the overall reaction as a function of pH (rows) and T (columns)
      A <- result$a12
      # Convert dimensionless affinity (A/2.303RT) to delta G (kJ / mol)
      # For temperature values to be applied correctly, we need to transpose
      # to get T into the rows of A (i.e., the first indexed dimension),
      # then transpose again to get back to the plot dimensions
      # For yvar == "P" this has no effect because TK is a scalar
      G.J <- t(convert(t(A), "G", T = TK))
      G.kJ <- G.J / 1000
      # Add color image
      image(a$vals$pH, a$vals[[yvar]], G.kJ, col = col, breaks = levels, add = TRUE)
      # Add contour lines
      contour(a$vals$pH, a$vals[[yvar]], G.kJ, levels = levels, labcex = 0.9, add = TRUE, method = method, lwd = lwd)

      # Replot tick marks
      thermo.axis()
      # Add legend
      legend(legend.x, legend.space, bg = "white")
      legend(legend.x, legend.expr, bty = "n")
      # Replot border
      box()
      # Add title
      # Use e.g. adenosine instead of adenosine_to_AMP
      short_reactant <- strsplit(reactant, "_")[[1]][1]
      # Change pyruvic acid to pyruvate
      short_reactant[short_reactant == "pyruvic acid"] <- "pyruvate"
      main <- bquote(log~italic(a)[.(short_reactant)]==.(loga_reactant))
      title(main, cex.main = 1.3)
      # Add panel label - outside x range
      label.figure(letters[iloga + dlab], cex = 1.6, xfrac = 0.03)

    }

  }

  # Put together a data frame for describe.reaction()
  # Use the first three basis species for the basic reaction: reactant, H3PO4, product
  coeff <- -as.numeric(species()[1:3])
  # For cAMP, double the reaction to consume one H3PO4
  if(reactant == "adenosine_to_cAMP") coeff <- coeff * 2
  formula <- names(species())[1:3]
  # For getting the name, use ispecies rather than formula (to avoid mixing up glucose and fructose)
  name <- info(basis()$ispecies[1:3])$name
  # Add H2O
  coeff <- c(coeff, 1)
  formula <- c(formula, "H2O")
  name <- c(name, "water")
  if(P_source != "P") {
    # For extended reactions, replace H3PO4 and H2O with P_source and P_remainder
    # Do names first
    name[formula == "H3PO4"] <- info(info(names(result$P_reaction)[1]))$name
    name[formula == "H2O"] <- info(info(names(result$P_reaction)[2]))$name
    # Then formulas
    formula[formula == "H3PO4"] <- names(result$P_reaction)[1]
    formula[formula == "H2O"] <- names(result$P_reaction)[2]
  }
  # Use ionized forms for names
  name[name == "H3PO4"] <- "Pi"
  name[name == "H4P2O7"] <- "PP"
  name[name == "acetylphosphate0"] <- "acetylphosphate"
  name[name == "acetic acid"] <- "acetate"
  name[name == "H2AMP"] <- "AMP"
  name[name == "H3ADP"] <- "ADP"
  name[name == "H4ATP"] <- "ATP"
  name[name == "H2UMP"] <- "UMP"
  name[name == "cyclic-HAMP"] <- "cyclic-AMP"
  name[name == "pyruvic acid"] <- "pyruvate"
  # Construct data frame
  reaction <- data.frame(coeff, formula, name)
  # Use names except for inorganic species
  use.name <- c(TRUE, TRUE, TRUE, TRUE)
  use.name[formula == "H2O"] <- FALSE
  # Uncomment these to use formulas instead of Pi and PP
  #use.name[formula == "H3PO4"] <- FALSE
  #use.name[formula == "H2P4O7"] <- FALSE
  iname <- which(use.name)
  # Change minus signs to short hyphens
  for(i in iname) reaction$name[i] <- hyphen.in.pdf(reaction$name[i])
  # Get expression for reaction 
  reaction_expr <- describe.reaction(reaction, iname = iname)
  # Add reaction to plot
  opar <- par(mar = c(0, 0, 0, 0))
  plot.new()
  text(0.5, 0.5, reaction_expr, cex = 1.4)
  par(opar)

  # After we're done with the plot, calculate standard transformed Gibbs energy (ΔG°') at pH 7
  result <- phosphorylate(reactant, P_source, const_pH = 7)
  # The affinity of the overall reaction
  A <- result$a12
  # Convert to Delta G
  TK <- convert(25, "K")
  G.J <- convert(A, "G", T = TK)
  G.kJ <- G.J / 1000
  print(paste("DeltaG0' at 25 degC and pH = 7 (kJ/mol):", round(G.kJ, 3)))
  # Return the calculated value
  G.kJ

} # end of phospho_plot()

