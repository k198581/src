source("lib/load_data_processing.R")
source("lib/load_scoring_matrix.R")
source("lib/load_exec_align.R")
source("lib/load_verif_lib.R")
source("psa/pairwise_pmi.R")
source("psa/pf-pmi.R")
source("verification/verification_psa.R")
source("parallel_config.R")

file <- "ansrate_pf-pmi"
dir <- "pairwise_pf-pmi"
ext = commandArgs(trailingOnly = TRUE)[1]
path <- MakePath(file, dir, ext)

# Create an itnitial scoring matrix and a list of PSAs.
list.words <- GetPathList()
s          <- MakeEditDistance(Inf)
psa.list   <- PSAforEachWord(list.words, s, dist = T)

# Update the scoring matrix using the PF-PMI.
pmi.o    <- PairwisePMI(psa.list, list.words, s, UpdatePFPMI)
pmi.mat  <- pmi.o$pmi.mat
s        <- pmi.o$s
psa.list <- pmi.o$psa.list

# Execute the PSA for each word.
VerificationPSA(psa.list, path$ansrate.file, path$output.dir)

# Save the matrix of the PMIs and the scoring matrix.
rdata.path <- MakeMatPath("matrix_psa_pf-pmi", "score_psa_pf-pmi", ext)
save(pmi.mat, file = rdata.path$rdata1)
save(s, file = rdata.path$rdata2)
