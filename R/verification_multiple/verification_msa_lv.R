source("lib/load_data_processing.R")
source("lib/load_scoring_matrix.R")
source("lib/load_exec_align.R")
source("lib/load_verif_lib.R")
source("verification_multiple/CalcAccMSA.R")
source("parallel_config.R")

ansrate <- "ansrate_msa_lv"
multiple <- "multiple_lv"
ext = commandArgs(trailingOnly=TRUE)[1]
path <- MakePath(ansrate, multiple, ext)

# Get the all of files path.
list.words <- GetPathList()
s <- MakeEditDistance(Inf)
msa.list <- MSAforEachWord(list.words, s)
CalcAccMSA(msa.list, list.words, path$ansrate.file, path$output.dir)
