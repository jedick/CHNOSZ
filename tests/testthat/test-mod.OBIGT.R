context("mod.OBIGT")

# clear out any previous basis definition or database alterations
suppressMessages(reset())

test_that("minimal usage of mod.OBIGT() creates usable data entries", {
  # we need at least a name and some property
  expect_error(mod.OBIGT("test"), "species name and a property")
  # a valid formula is needed
  expect_warning(expect_error(mod.OBIGT("test", date=today()), "is not a simple chemical formula"),
               "please supply a valid chemical formula")
  # the default state is aq
  expect_message(itest <- mod.OBIGT("test", formula="Z0", date=today()), "added test\\(aq\\)")
  # set the charge so following test use hkf() rather than AkDi()
  mod.OBIGT("test", z = 0)
  # we should get NA values of G for a species with NA properties 
  expect_true(all(is.na(subcrt(itest)$out[[1]]$G)))
  # a single value of G comes through to subcrt
  mod.OBIGT("test", G=100)
  expect_equal(subcrt("test", T=25, P=1)$out[[1]]$G, 100)
  # values for Cp and c1 integrate to the same values of G
  G.Cp <- subcrt(mod.OBIGT(list(name="test", S=0, Cp=100)))$out[[1]]$G
  G.c1 <- subcrt(mod.OBIGT(list(name="test", S=0, c1=100)))$out[[1]]$G
  expect_equal(G.Cp, G.c1)
})

test_that("mod.OBIGT() works with numeric argument (species index)", {
  # test added 20200707
  i1 <- info("CH4")
  # previously this failed with an error
  i2 <- mod.OBIGT(i1, name = "methane")
  expect_identical(i1, i2)
  expect_identical(info(i2)$name, "methane")
})
