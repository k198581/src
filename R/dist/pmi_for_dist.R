source("parallel_config.R")

source("dist/tools/get_regions.R")

source("lib/load_scoring_matrix.R")
source("lib/load_nwunsch.R")

source("psa/pmi.R")
source("psa/pairwise_pmi.R")

load("all_list.RData")

LD <- function(c1, c2, s) {
  
  N1 <- length(c1)
  N2 <- length(c2)
  
  score.vec <- NULL
  psa.list  <- list()
  k <- 1
  for (i in 1:N1) {
    for (j in 1:N2) {
      psa    <- NeedlemanWunsch(c1[[i]], c2[[j]], s, select.min = T)
      as     <- psa$score
      c1.len <- length(c1[[i]]) - 1
      c2.len <- length(c2[[j]]) - 1
      ldn    <- as / max(c1.len, c2.len)
      
      psa$score     <- ldn
      psa.list[[k]] <- psa
      score.vec     <- c(score.vec, ldn)
      
      k <- k + 1
    }
  }
  psa.list[which(score.vec == min(score.vec))[1]]
}

PSAforEachConcept <- function(r1, r2, s) {
  
  n1 <- NULL
  n2 <- NULL
  
  N <- length(r1)
  for (i in 1:N) {
    n1[i] <- length(r1[[i]])
    n2[i] <- length(r2[[i]])
  }
  
  n <- n1 * n2
  zero <- rev(which(n == 0))
  
  for (i in zero) {
    r1 <- r1[-i]
    r2 <- r2[-i]
  }
  N <- length(r1)
  
  score.vec <- NULL
  concepts  <- names(r1)
  psa.list  <- list()
  k <- 1
  for (i in 1:N) {
    for (j in 1:N) {
      psa           <- unlist(LD(r1[[i]], r2[[j]], s), recursive = F)
      psa.list[[k]] <- psa
      score.vec     <- c(score.vec, psa$score)
      
      k <- k + 1
    }
  }
  psa.list[which(score.vec == min(score.vec))[1]]
}

PSAforEachResion <- function(all.list, s) {
  
  # Make the PSA list.
  r <- t(combn(95, 2))
  N <- dim(r)[1]
  
  psa.list <- foreach (i = 1:N) %dopar% {
    
    k <- r[i, 1]
    l <- r[i, 2]
    
    # region
    r1 <- all.list[[k]]
    r2 <- all.list[[l]]
    
    PSAforEachConcept(r1, r2, s)
  }
  
  psa.list
}

# Find for the PSA which is lowest score.
find_the_lowest <- function(psa.list) {
  M <- length(psa.list)
  for (i in 1:M) {
    N <- length(psa.list[[i]])
    vec <- foreach (j = 1:N, .combine = "c") %dopar% {
      psa.list[[i]][[j]]$score
    }
    psa.list[[i]] <- psa.list[[i]][which(vec == min(vec))]
  }
  psa.list
}

# Update the scoring matrix using the PMI.
s <- MakeEditDistance(Inf)  # the initial scoring matrix
cat("\n")
print("Initial PSA")
psa.list <- PSAforEachResion(all.list, s)  # the initial alignments
#print("Find the lowest")
#psa.list <- find_the_lowest(psa.list) 
cat("\n")

s.old <- s
N <- length(s.old)
for (i in 1:N) {
  s.old[i] <- 0
}
# START OF LOOP
i <- 0
while(1) {
  i <- i + 1
  cat("\n")
  print(paste("loop:", i))
  diff <- N - sum(s == s.old)
  if (diff == 0) break
  # Compute the new scoring matrix that is updated by the PMI-weighting.
  print("Updating the PMI.")
  s.old <- s
  rlt.pmi <- UpdatePMI(psa.list, s)
  pmi.mat <- rlt.pmi$pmi.mat
  s <- rlt.pmi$s
  # Compute the new PSA using the new scoring matrix.
  print("PSAforEachResion")
  psa.list <- PSAforEachResion(all.list, s)
  #print("find_the_lowest")
  #psa.list <- find_the_lowest(psa.list) 
}

save(pmi.mat, file = "pmi_mat.RData")
save(s, file = "pmi_score.RData")

print("Finished!!")
