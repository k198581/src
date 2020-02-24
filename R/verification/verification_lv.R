source("verification/VerificationPSA.R")
source("lib/load_scoring_matrix.R")
source("parallel_config.R")

VerificationLV <- function(file, dir, ext=NULL) {
  
  # Set the path of the matching rate.
  ansrate.file <- paste("../../Alignment/", file, "_", ext, "_", format(Sys.Date()), ".txt", sep = "")
  
  # Set the path of the PSA directory.
  output.dir <- paste("../../Alignment/", dir, "_", ext, "_", format(Sys.Date()), "/", sep = "")
  if (!dir.exists(output.dir))
    dir.create(output.dir)
  
  # Make the scoring matrix.
  s <- MakeEditDistance(Inf)
  
  # Execute the PSA for each word.
  VerificationPSA(ansrate.file, output.dir, s)
  
}

file <- "ansrate_lv"
dir <- "pairwise_lv"
ext = commandArgs(trailingOnly=TRUE)[1]
VerificationLV(file, dir, ext)
