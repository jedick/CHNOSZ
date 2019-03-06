context("retrieve")

test_that("packaged stoichiometric matrix matches the default database", {
  expect_identical(rownames(thermo()$stoich), thermo()$obigt$formula)
  # test added 20190303
  # if this test fails, update extdata/thermo/stoich.csv.xz and rebuild the package:
  # reset()
  # formula <- thermo()$obigt$formula
  # stoich <- i2A(formula)
  # write.csv(stoich, "stoich.csv")
  # system("xz -f stoich.csv")
})

test_that("errors and recalculations produce expected messages", {
  expect_error(retrieve(c("A", "B", "C")), '"A" is not an element')
  expect_error(retrieve(c("A", "B", "C", "D")), '"A", "D" are not elements')
  add.obigt("SUPCRT92")
  expect_message(retrieve("Ti"), "updating stoichiometric matrix")
  reset()
})

test_that("retrieve('all') works", {
  all1 <- sort(retrieve("all"))
  all2 <- sort(retrieve(as.list(colnames(thermo()$stoich))))
  expect_identical(all1, all2)
})
