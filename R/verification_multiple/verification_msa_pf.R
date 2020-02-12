source("lib/load_data_processing.R")
source("lib/load_scoring_matrix.R")
source("lib/load_exec_align.R")
source("verification_multiple/CalcAccMSA.R")
source("parallel_config.R")

ansrate <- "ansrate_msa_pf"
multiple <- "multiple_pf"

list.words <- GetPathList()  # the list of words
s <- MakeEditDistance(Inf)  # the scoring matrix

for (pen in (-1)) {
  
  # matchingrate path
  ansrate.file <- paste("../../Alignment/", ansrate, "_", pen, "_", format(Sys.Date()), ".txt", sep = "")
  
  # result path
  output.dir <- paste("../../Alignment/", multiple, "_", pen, "_", format(Sys.Date()), "/", sep = "")
  if (!dir.exists(output.dir)) {
    dir.create(output.dir)
  }

  # Make the list of the MSAs.
  msa.list <- MSAforEachWord(list.words, s)
  # Calculate the accuracy of the MSAs.
  CalcAccMSA(msa.list, list.words, ansrate.file, output.dir)
  
}
