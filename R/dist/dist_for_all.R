source("parallel_config.R")

source("dist/tools/get_regions.R")
source("dist/tools/make_mat.R")
source("dist/tools/dist.R")

load("all_list.RData")
load("i.vec.RData")
load("j.vec.RData")

dist_for_all <- function(method="lv") {
  
  r.pairs <- t(combn(95, 2))
  N <- dim(r.pairs)[1]
  
  psa.list <- foreach (i = j.vec) %dopar% {
    
    k <- r.pairs[i, 1]
    l <- r.pairs[i, 2]
    
    # region
    r1 <- all.list[[k]]
    r2 <- all.list[[l]]
    
    mat.o  <- MakeMat(r1, r2, method)
    
    dist        <- list()
    dist$method <- method
    dist$pair   <- paste(regions[k], regions[l])  # for using UTF-8.
    dist$pair   <- unlist(strsplit(dist$pair, " "))
    
    dist$psa.list <- mat.o$psa.list
    return(dist)
  }
  
  psa.list
}

file   = commandArgs(trailingOnly=TRUE)[1]
method = commandArgs(trailingOnly=TRUE)[2]

if (is.na(method)) {
  method <- "lv"
} else {
 # NOP 
}

if (is.na(file)) {
  stop("An output file name must be selected.")
} else {
  print(paste("File:", file))
  print(paste("Method:", method))
  psa.list <- dist_for_all(method)
  save(psa.list, file = paste("psa_list_", file, ".RData", sep = ""))
}

print("Finished!!")
