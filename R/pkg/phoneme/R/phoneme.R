# Consonants
C <- as.vector(read.table("R/symbols/consonants.txt")[, 1])

# Vowels
V <- as.vector(read.table("R/symbols/vowels.txt")[, 1])

# Consonant features
mat.C.feat <- as.matrix(read.table("R/features/consonants.txt"))
dimnames(mat.C.feat) <- list(C, NULL)

# Vowel features
mat.V.feat <- as.matrix(read.table("R/features/vowels.txt"))
dimnames(mat.V.feat) <- list(V, NULL)

for (j in 1:5) {
  mat.C.feat[, j] <- paste("C", mat.C.feat[, j], j, sep = "")
  mat.V.feat[, j] <- paste("V", mat.V.feat[, j], j, sep = "")
}

C.feat <- unique(as.vector(mat.C.feat))
V.feat <- unique(as.vector(mat.V.feat))

CV <- c(C, V)
CV.feat <- c(C.feat, V.feat)
mat.CV.feat <- rbind(mat.C.feat, mat.V.feat)

if (0) {
  usethis::use_data(C, C.feat, mat.C.feat,
                    V, V.feat, mat.V.feat,
                    CV, CV.feat, mat.CV.feat, internal = T)
}
