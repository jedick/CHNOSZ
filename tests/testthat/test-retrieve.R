context("retrieve")

test_that("packaged stoichiometric matrix matches the default database", {
  expect_identical(rownames(thermo()$stoich), thermo()$obigt$formula)
  # test added 20190303
  # if this test fails, update extdata/thermo/stoich.csv.xz and rebuild the package:
  # reset()
  # formula <- thermo()$obigt$formula
  # stoich <- i2A(formula)
  # write.csv(stoich, "stoich.csv")
  # system("xz stoich.csv")
})

test_that("errors and recalculations produce expected messages", {
  # this should give an error about one non-element
  expect_error(retrieve(c("A", "B", "C")), '"A" is not an element')
  # this should give an error about two non-elements
  expect_error(retrieve(c("A", "B", "C", "D")), '"A", "D" are not elements')
  # this should recalculate the stoichiometric matrix
  add.obigt("SUPCRT92")
  expect_message(retrieve("Ti"), "creating stoichiometric matrix")
  reset()
})
