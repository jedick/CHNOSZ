# Load default settings for CHNOSZ
reset()

info <- "expr.species() produces expected errors"
expect_error(expr.species(c("H2O", "CO2")), "more than one species", info = info)

# 20220608
info <- "expr.species() handles non-integer coefficients"
es <- expr.species("FeS1.33")
expect_equal(es[[length(es)]], "1.33", info = info)

info <- "expr.species() handles non-integer charge"
es <- expr.species("FeS+1.33")
expect_equal(es[[length(es)]], "+1.33", info = info)

info <- "describe.property() does not accept NULL values"
expect_error(describe.property(), "property or value is NULL", info = info)

# Tests added on 20260109

info <- "describe.basis() runs without error"
basis("CHNOS+")
basis("FeCHNOS+")
basis("O2", "HM")
expect_silent(describe.basis(), info = info)
expect_silent(describe.basis(oneline = TRUE), info = info)

info <- "describe.property() runs without error"
property <- c("P", "T", "Eh", "pH", "IS")
value <- c(1, 42.42, -1, 7, 0.1)
expect_silent(describe.property(property, value), info = info)
expect_silent(describe.property(property, value, oneline = TRUE), info = info)

info <- "ratlab() and syslab() run without error"
expect_silent(ratlab(), info = info)
expect_silent(syslab(), info = info)

info <- "split.formula() works with H-citrate-2"
expect_silent(result <- CHNOSZ:::split.formula("H-citrate-2"), info = info)
expect_equivalent(result["H-citrate"], 1, info = info)
expect_equivalent(result["Z"], -2, info = info)

info <- "split.formula() works with 2-octanone"
expect_equal(CHNOSZ:::split.formula("2-octanone"), "2-octanone")
