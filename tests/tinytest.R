# Changed to tinytest 20220129

if ( requireNamespace("tinytest", quietly = TRUE) ){
  # Use 4-number package versioning while developing and 3-number versioning for CRAN release
  # https://stackoverflow.com/questions/36166288/skip-tests-on-cran-but-run-locally
  # (mentioned in tinytest/vignettes/using_tinytest.pdf)
  at_home <- length(unclass(packageVersion("CHNOSZ"))[[1]]) == 4
  tinytest::test_package("CHNOSZ", at_home = at_home)
}
