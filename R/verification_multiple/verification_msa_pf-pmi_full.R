source("lib/load_data_processing.R")
source("lib/load_verif_lib.R")
source("lib/load_scoring_matrix.R")
source("lib/load_exec_align.R")
source("msa/ProgressiveAlignment.R")
source("msa/BestFirst.R")
source("psa/pairwise_pmi.R")
source("psa/pairwise_pf-pmi.R")
source("verification_multiple/change_list_msa2psa.R")
source("verification_multiple/VerificationMSA.R")
source("parallel_config.R")

ansrate <- "ansrate_msa"
multiple <- "multiple"

# matchingrate path
ansrate.file <- paste("../../Alignment/", ansrate, "_pf-pmi_", format(Sys.Date()), ".txt", sep = "")

# result path
output.dir <- paste("../../Alignment/", multiple, "_pf-pmi_", format(Sys.Date()), "/", sep = "")
if (!dir.exists(output.dir)) {
  dir.create(output.dir)
}

# Compute the scoring matrix using the PMI method.
list.words <- GetPathList()
s <- MakeEditDistance(Inf)
psa.list <- PSAforEachWord(list.words, s)
s <- PairwisePFPMI(psa.list, list.words, s)
#save(s, file="scoring_matrix_msa_pmi.RData")

# For progressive
s.old <- s
N <- length(s.old)
for (i in 1:N) {
  s.old[i] <- 0
}

msa.list <- list()
while (1) {
  diff <- N - sum(s == s.old)
  if (diff == 0) break
  #
  for (w in list.words) {
    # Make the word list.
    gold.list <- MakeWordList(w["input"])  # gold alignment
    seq.list <- MakeInputSeq(gold.list)  # input sequences
    
    # Computes the MSA using the BestFirst method.
    print(paste("Start:", w["name"]))
    id <- as.numeric(w["id"])
    msa.init <- ProgressiveAlignment(psa.list[[id]], seq.list, s, F)
    msa.list[[id]] <- list()
    msa.list[[id]] <- BestFirst(msa.init, s, F)
    print(paste("End:", w["name"]))
  } 
  psa.list <- ChangeListMSA2PSA(msa.list, s)
  s.old <- s
  s <- PairwisePFPMI(psa.list, list.words, s)
}

# For verification
CalcAccMSA(msa.list, list.words, ansrate.file, output.dir)
