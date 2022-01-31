# Load default settings for CHNOSZ
reset()

# Test added 20190303
# If this test fails, update extdata/thermo/stoich.csv.xz and rebuild the package:
# reset()
# formula <- thermo()$OBIGT$formula
# stoich <- i2A(formula)
# write.csv(stoich, "stoich.csv")
# system("xz -f stoich.csv")
info <- "Prebuilt stoichiometric matrix is consistent with the default database"
expect_identical(rownames(thermo()$stoich), thermo()$OBIGT$formula, info = info)

info <- "Errors and recalculations produce expected messages"
expect_error(retrieve(c("A", "B", "C")), '"A" is not an element', info = info)
expect_error(retrieve(c("A", "B", "C", "D")), '"A", "D" are not elements', info = info)
add.OBIGT("SUPCRT92")
expect_message(retrieve("Ti"), "updating stoichiometric matrix", info = info)

info <- "retrieve('all') works"
all1 <- sort(retrieve("all"))
all2 <- sort(retrieve(as.list(colnames(thermo()$stoich))))
expect_identical(all1, all2, info = info)
