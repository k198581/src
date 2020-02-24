source("lib/load_data_processing.R")
source("lib/load_scoring_matrix.R")
source("lib/load_exec_align.R")
source("verification_multiple/CalcAccMSA.R")
source("parallel_config.R")

VerificationMSALV <- function(ansrate, multiple, ext=NULL) {

  # matchingrate path
  ansrate.file <- paste("../../Alignment/", ansrate, "_", ext, "_", format(Sys.Date()), ".txt", sep = "")
  
  # result path
  output.dir <- paste("../../Alignment/", multiple, "_", ext, "_", format(Sys.Date()), "/", sep = "")
  if (!dir.exists(output.dir)) {
    dir.create(output.dir)
  }
  
  # Get the all of files path.
  list.words <- GetPathList()
  s <- MakeEditDistance(Inf)
  msa.list <- MSAforEachWord(list.words, s)
  CalcAccMSA(msa.list, list.words, ansrate.file, output.dir)
  
}

ansrate <- "ansrate_msa_lv"
multiple <- "multiple_lv"
ext = commandArgs(trailingOnly=TRUE)[1]
VerificationMSA(ansrate, multiple, ext)
