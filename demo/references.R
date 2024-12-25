# CHNOSZ/demo/reference.R
## Cross-checking references
library(CHNOSZ)

# The available reference keys
available_refs <- thermo()$refs$key

# The used reference keys
# List all files in OBIGT and OBIGT/testing
OBIGT_files <- dir(system.file("extdata/OBIGT", package = "CHNOSZ"), pattern = ".csv", full.names = TRUE)
testing_files <- dir(system.file("extdata/OBIGT/testing", package = "CHNOSZ"), pattern = ".csv", full.names = TRUE)
all_files <- c(OBIGT_files, testing_files)
all_data <- lapply(all_files, read.csv)
all_data <- do.call(rbind, all_data)
used_refs <- na.omit(c(all_data$ref1, all_data$ref2))

# Print messages for unavailable and unused refs
unavailable_refs <- used_refs[!used_refs %in% available_refs]
unused_refs <- available_refs[!available_refs %in% used_refs]
if(length(unavailable_refs) > 0 | length(unused_refs) > 0) {

  if(length(unavailable_refs) > 0) {
    print("References used in database but not available in refs.csv:")
    print(unique(unavailable_refs))
  }

  if(length(unused_refs) > 0) {
    print("References available in refs.csv but not used in database:")
    print(unique(unused_refs))
  }

} else {
  
  print("All references accounted for!")

}
