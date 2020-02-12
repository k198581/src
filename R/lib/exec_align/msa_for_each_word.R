source("lib/load_data_processing.R")
source("lib/load_verif_lib.R")
source("msa/ProgressiveAlignment.R")
source("msa/BestFirst.R")

MSAforEachWord <- function(list.words, s, similarity=F) {
  # Compute the MSA for each word.
  # Args:
  #   ansrate.file: The path of the matching rate file.
  #   output.dir:   The path of the MSA directory.
  #   s:            The scoring matrix.
  #
  # Returns:
  #   The list of MSA for each word.
  
  msa.list <- list()
  for (w in list.words) {
    
    # Make the word list.
    gold.list <- MakeWordList(w["input"])  # gold alignment
    seq.list <- MakeInputSeq(gold.list)     # input sequences
    
    # Computes the MSA using the BestFirst method.
    print(paste("Start:", w["name"]))
    psa.list <- MakePairwise(seq.list, s, select.min=!similarity)
    psa.init <- ProgressiveAlignment(psa.list, seq.list, s, similarity)
    msa.list[[as.numeric(w["id"])]] <- list()
    msa.list[[as.numeric(w["id"])]] <- BestFirst(psa.init, s, similarity)
    print(paste("End:", w["name"]))
    
  }
  
  return(msa.list)
}
