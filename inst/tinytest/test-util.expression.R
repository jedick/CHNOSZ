# Load default settings for CHNOSZ
reset()

info <- "expr.species() produces expected errors"
expect_error(expr.species(c("H2O", "CO2")), "more than one species", info = info)

# 20220608
info <- "expr.species() handles non-integer coefficients"
es <- expr.species("FeS1.33")
expect_equal(es[[length(es)]], "1.33")

info <- "expr.species() handles non-integer charge"
es <- expr.species("FeS+1.33")
expect_equal(es[[length(es)]], "+1.33")
