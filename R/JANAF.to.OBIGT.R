# JANAF.to.OBIGT.R
# Convert JANAF data to OBIGT format
# 20260604 jmd

# The first four lines of https://janaf.nist.gov/tables/C-127.txt
# Ethyne (C2H2)   C2H2(g)
# T(K)    Cp      S       -[G-H(Tr)]/T    H-H(Tr) delta-f H       delta-f G       log Kf
# 0       0.000   0.000   INFINITE        -10.012 227.288 227.288 INFINITE
# 100     29.347  163.294 234.338 -7.104  227.078 221.011 -115.444

JANAF.to.OBIGT <- function(file, abbrv = NULL, T_max = 1500, MAE_max = 1, plot_Cp = FALSE) {

  # Preprocess file: replace strings to get equal numbers of columns for read.table()
  lines <- readLines(file)
  lines <- gsub("Cp LAMBDA MAXIMUM", "NA\tNA\tNA", lines)
  lines <- gsub("TRANSITION", "NA\tNA\tNA", lines)
  tmpfile <- tempfile(fileext=".txt")
  writeLines(lines, tmpfile)

  # Process abbrv argument
  abbrv = ifelse(is.null(abbrv), NA_character_, abbrv)
  # Read the table, skipping the first three lines
  dat <- read.table(tmpfile, skip = 3)
  # Add column names
  colnames(dat) <- c("T", "Cp", "S", "DG", "DH", "H", "G", "logK")
  # Read the first line to get name, formula, and state
  first <- readLines(tmpfile, n = 1)
  name <- tolower(strsplit(first, " ")[[1]][1])
  # Get formula with state, e.g. C2H2(g)
  formula_state <- strsplit(first, "\t")[[1]][2]
  # Get formula by itself, e.g. C2H2
  formula <- strsplit(formula_state, "\\(")[[1]][1]
  # Get state (inside parentheses)
  m <- regexec("\\(([^)]*)\\)", formula_state)
  state_orig <- regmatches(formula_state, m)[[1]][2]
  # Convert g to gas
  # TODO: work with other states (aq, cr)
  state <- switch(state_orig, g = "gas", NA)
  if(is.na(state)) stop(paste("unrecognized state:", state_orig))
  # Format today's date for ISO 8601, e.g. 2026-06-04
  date <- format(Sys.time(), "%Y-%m-%d")
  # Find row for 298.15 K
  i298 <- dat$T==298.15

  # Find rows to use for Cp values (100 to 1500 K)
  iCp <- dat$T >= 100 & dat$T <= T_max
  Cp <- dat$Cp[iCp]
  # Use linear model to fit Cp equation
  # Cp = a + b*T + c*T^-2 + d*T^-0.5 + e*T^2
  T <- dat$T[iCp]
  T_2 <- T ^ -2
  T_0.5 <- T ^ -0.5
  T2 <- T ^ 2
  Cp_lm <- lm(Cp ~ T + T_2 + T_0.5 + T2)
  Cp_predicted <- predict(Cp_lm)

  if(plot_Cp) {
    # Plot actual and predicted values
    plot(T, Cp, xlab = quote(italic(T)~"(K)"), ylab = axis.label("Cp"))
    lines(T, Cp_predicted)
    legend("topleft", c("actual", "predicted"), pch = c(1, NA), lty = c(0, 1))
    title(formula_state, font.main = 1)
  }
  # Calculate MAE
  MAE <- mean(abs(Cp_predicted - Cp))
  print(paste0("MAE for Cp of ", name, " ", formula_state, ": ", round(MAE, 2)))
  # Error if MAE > max MAE
  stopifnot(MAE <= MAE_max)

  # Put together parameters
  PAR <- list(
    name = name,
    abbrv = abbrv,
    formula = formula,
    state = state,
    ref1 = "JANAF98",
    ref2 = "OBIGT26",
    date = date,
    model = "CGL",
    E_units = "J",
    G = dat$G[i298] * 1000,
    H = dat$H[i298] * 1000,
    S = dat$S[i298],
    Cp = dat$Cp[i298],
    V = 0,
    a = signif(Cp_lm$coefficients[1], 5),
    b = signif(Cp_lm$coefficients[2], 5),
    c = signif(Cp_lm$coefficients[3], 5),
    d = signif(Cp_lm$coefficients[4], 5),
    e = signif(Cp_lm$coefficients[5], 5),
    f = 0,
    lambda = 0,
    T = 1500,
    zap = TRUE
  )

  # Add data to OBIGT
  inew <- do.call(mod.OBIGT, PAR)
  return(inew)

}
