library(foreach)
library(doParallel)
registerDoParallel(detectCores())

source("lib/load_data_processing.R")
#source("verification/VerificationPSA.R")
source("lib/load_verif_lib.R")
source("lib/load_scoring_matrix.R")
source("psa/pmi.R")
source("psa/pf-pmi.R")

# get the all of files path
filesPath <- GetPathList()

ansrate <- "ansrate_pf-pmi"
pairwise <- "pairwise_pf-pmi"

files <- GetPathList()
input.list <- MakeInputList(files)

load("../../Alignment/scoring_matrices.RData")

# matchingrate path
ansrate.file <- paste("../../Alignment/", ansrate, "_", format(Sys.Date()), ".txt", sep = "")

# result path
output.dir <- paste("../../Alignment/", pairwise, "_", format(Sys.Date()), "/", sep = "")
if (!dir.exists(output.dir))
  dir.create(output.dir)

# conduct the alignment for each files
foreach.rlt <- foreach (f = filesPath) %dopar% {
  
  # make the word list
  gold.list <- MakeWordList(f["input"])
  input.list <- MakeInputSeq(gold.list)
  
  # making the gold standard alignments
  gold.aln <- MakeGoldStandard(gold.list)

  psa.aln <- list()
  for (p in 1:5)
    psa.aln[[p]] <- MakePairwise(input.list, s.list[[p]], select.min = T)

  n <- length(psa.aln[[1]])
  score <- c(0, 0, 0, 0, 0)
  for (p in 1:5)
    for (i in 1:n)
      score[p] <- score[p] + psa.aln[[p]][[i]]$score
  
  score.min <- which(score==min(score))[1]
  
  psa.aln <- psa.aln[[score.min]]
  
  #######
  # calculating the matching rate
  matching.rate <- VerifAcc(psa.aln, gold.aln)
  # output gold standard
  OutputAlignment(f["name"], output.dir, ".lg", gold.aln)
  # output pairwise
  OutputAlignment(f["name"], output.dir, ".aln", psa.aln)
  # output match or mismatch
  OutputAlignmentCheck(f["name"], output.dir, ".check", psa.aln, gold.aln)
  
  # Returns the matching rate to the list of foreach.
  c(f["name"], matching.rate)
}

# Outputs the matching rate
matching.rate.mat <- list2mat(foreach.rlt)
matching.rate.mat <- matching.rate.mat[order(matching.rate.mat[, 1]), , drop=F]
write.table(matching.rate.mat, ansrate.file, quote = F)
