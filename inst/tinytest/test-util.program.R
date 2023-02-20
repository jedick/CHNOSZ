# This is a long test ... only run it "at home" 20220131
if(!at_home()) exit_file("Skipping long test")

# Load default settings for CHNOSZ
reset()

# These tests show inefficient uses of parallelization
# (overhead is greater than the savings from multiple cores).
# They are here just to test that the functions are working.

info <- "palply() launches calculations on multiple cores"
if(min(getOption("mc.cores"), 2) > 1 & parallel::detectCores() > 1) {
  x <- 1:1001
  # for this we don't have to export any variables so varlist == ""
  expect_message(palply("", 1:length(x), function(i) i^2), "running 1001 calculations", info = info)
}

info <- "Other functions are calling palply() properly"
if(min(getOption("mc.cores"), 2) > 1 & parallel::detectCores() > 1) {
  basis("CHNOS")
  ip <- 1:nrow(thermo()$protein)
  expect_message(a <- affinity(iprotein = rep(ip, 3)), "affinity running .* calculations", info = info)
  expect_message(e <- equilibrate(a, normalize = TRUE), "equil.boltzmann running .* calculations", info = info)
  # Test reaction method
  species(c("CO2", "acetic acid"))
  a <- affinity(O2 = c(-90, -60, 1000))
  expect_message(e <- equilibrate(a), "equil.reaction running 1000 calculations", info = info)
}
