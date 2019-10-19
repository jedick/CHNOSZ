context("wjd")

test_that("wjd() gives results similar to White et al., 1958", {
  # the values from last column of Table III in the paper
  X <- c(0.040668, 0.147730, 0.783153, 0.001414, 0.485247, 0.000693, 0.027399, 0.017947, 0.037314, 0.096872)
  w <- wjd()
  expect_equal(X, w$X, tolerance=1e-4)
})

test_that("guess() operates on intermediate compositions but fails with endmembers", {
  alkanes <- c("hexane", "heptane", "octane", "nonane")
  ialk <- info(alkanes, "liq")
  expect_error(run.guess(ialk, "C6H14"), "there are only 0")
  # hmm, on windows this has a length of 5 (20120626)
  # probably should filter out guesses with very low abundances
  #expect_true(length(run.guess(ialk, "C7H16"))==4)
  expect_true(length(run.guess(ialk, "C8H18"))==5)
  expect_error(run.guess(ialk, "C9H20"), "there are only 0")
})

# references

# White, W. B., Johnson, S. M. and Dantzig, G. B. (1958) 
#   Chemical equilibrium in complex mixtures. 
#   J. Chem. Phys. 28, 751--755. https://doi.org/10.1063/1.1744264
