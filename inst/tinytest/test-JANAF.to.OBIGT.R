# Load default settings for CHNOSZ
reset()

# Tests added on 20260605

info <- "Error message for unrecognized state"
file <- system.file("extdata/misc/Fe-001.txt", package = "CHNOSZ")
expect_error(inew <- JANAF.to.OBIGT(file), "unrecognized state: cr", info = info)

info <- "Works as expected for ethyne(gas)"
# We need to set the current date in OBIGT so no changes are detected below
mod.OBIGT("ethyne", state = "gas", date = format(Sys.time(), "%Y-%m-%d"))
file <- system.file("extdata/misc/C-127.txt", package = "CHNOSZ")
expect_message(inew <- JANAF.to.OBIGT(file, abbrv = "acetylene"), "no change for ethyne\\(gas\\)", info = info)

info <- "Error message for too high MAE"
expect_error(inew <- JANAF.to.OBIGT(file, abbrv = "acetylene", MAE_max = 0.1), "MAE <= MAE_max is not TRUE")
