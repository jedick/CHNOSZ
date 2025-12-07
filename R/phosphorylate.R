# Calculate affinity of phosphorylation reactions taking account of speciation
# 20251206 first version (extracted from sugars paper script) jmd
phosphorylate <- function(reactant, P_source, loga_reactant = 0, loga_product = 0, loga_P_source = 0, loga_P_remainder = 0, const_pH = 7, ...) {

  # Affinity is calculated by summing:
  # 1) reactant + P = product + H2O (m1) 
  # 2) P_source + H2O = P_remainder + P (m2)
  # NOTE: reaction 2 is defined in reverse (to form H2O), then the opposite of the affinity is used for the sum

  # Terminology:
  # Complex species are species that take multiple forms (e.g. ionization or oxidation states) depending on pH or other variables
  # A basic reaction has complex species that are valid basis species (independent components)
  # --> This is what mosaic() deals with
  # An extended reaction has complex species all of which cannot be used as basis species (they are non-independent)
  # However, the energetics of an extended reaction can be modeled as a sum of basic reactions
  # --> This script has a function for this called extended_mosaic()

  # Examples:
  # Basic reaction:     acetic acid + P = acetylphosphate + H2O
  # Extended reaction:  acetic acid + PP = acetylphosphate + P
  # Complex species:    acetic acid = [acetic acid + acetate]

  ## Setup reaction 1
  if(reactant == "acetic acid") {
    # Basic reaction: acetic acid + P = acetylphosphate + H2O
    # Load initial species for mosaic reaction (uncharged species)
    basis(c("acetic acid", "H3PO4", "acetylphosphate0", "O2", "H+"))
    # The basis species we will swap through for mosaic
    bases <- list(
      c("H3PO4", "H2PO4-", "HPO4-2", "PO4-3"),
      c("acetylphosphate0", "acetylphosphate-1", "acetylphosphate-2", "acetylphosphate-3"),
      c("acetic acid", "acetate")
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
      c("H2AMP", "HAMP-", "AMP-2"),
      c("H3ADP", "H2ADP-", "HADP-2", "ADP-3")
    )
  } else if(reactant == "ADP") {
    # Basic reaction: ADP + P = ATP + H2O
    basis(c("H3ADP", "H3PO4", "H4ATP", "N2", "O2", "H+"))
    bases <- list(
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
      c("phosphoenolpyruvate", "phosphoenolpyruvate-1", "phosphoenolpyruvate-2", "phosphoenolpyruvate-3"),
      c("pyruvic acid", "pyruvate")
    )
  } else {
    stop(paste("unrecognized reactant:", reactant))
  }

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
  # Don't use species("acetylphosphate0"), because that reaction is just acetylphosphate0 = acetylphosphate0
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

}
