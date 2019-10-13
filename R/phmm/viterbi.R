Viterbi <- function(O, params, E) {
  # Computes the likelihood using forward algorithm.
  # 
  # Args:
  #   O: The list of the observation sequences and the lengths.
  #   params: The list of the parameters.
  #
  # Returns:
  #   The list of the distance matrices in each state.
  O1 <- O$O1
  O2 <- O$O2
  U <- O$U
  V <- O$V
  
  delta <- params["delta"]
  tau.M <- params["tau.M"]
  epsilon <- params["epsilon"]
  lambda <- params["lambda"]
  tau.XY <- params["tau.XY"]
  
  p.xy <- E$M
  q.x <- E$X
  q.y <- E$Y
  
  m.m <- 1-2*delta-tau.M
  xy.m <- 1-epsilon-lambda-tau.XY
  
  # Initialization
  v.M <- matrix(0, U+1, V+1)
  v.X <- matrix(0, U+1, V+1)
  v.Y <- matrix(0, U+1, V+1)
  
  v.M[2, 2] <- m.m
  v.X[2, 2] <- delta
  v.Y[2, 2] <- delta
  
  # Induction
  for (i in 2:U) {
    for (j in 2:V) {
      if ((i != 2) && (j != 2)) {
        v.M[i, j] <- p.xy[O1[i], O2[j]] * max(m.m * v.M[i-1, j-1], xy.m * v.X[i-1, j-1], xy.m * v.Y[i-1, j-1])
        v.X[i, j] <- q.x[1, O1[i]] * max(delta * v.M[i-1, j], epsilon * v.X[i-1, j], lambda * v.Y[i-1, j])
        v.Y[i, j] <- q.y[1, O2[j]] * max(delta * v.M[i, j-1], epsilon * v.Y[i, j-1], lambda * v.X[i, j-1])
      }
    }
  }
  
  Q <- c(NULL)
  args <- character()
  for (i in 2:U) {
    for (j in 2:V) {
      args <- c(v.M[i, j], v.X[i, j], v.Y[i, j])
      Q <- c(Q, head(which(args == max(args)), 1))
    }
  }
  
  Q[Q == 1] <- "M"
  Q[Q == 2] <- "X"
  Q[Q == 3] <- "Y"
  
  # Termination
  # Po <- tau.M * v.M[n, m] + tau.XY * (v.X[n, m] + v.Y[n, m])
  
  #v.val <- list(v.M, v.X, v.Y)
  #names(v.val) <- c(S)
  
  return(Q)
}
