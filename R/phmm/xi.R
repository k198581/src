Xi <- function(u, v, i, j, O, fb, A, E) {
  # Computes the transition probability from state i to j at u and j.
  #
  # Args:
  #   u, v: The values of positions on both of the observation sequences.
  #   i, j: The values of states.
  #   O: The list of the observation sequences and the lengths.
  #   fb: The list of the forward-backward algorithm.
  #   A: The matrix of transition probability.
  #   E: The list of the emission probability matrices.
  #
  # Returns:
  #   The value of the transition probability from state i to j at u and j.
  O1 <- O$O1
  O2 <- O$O2
  
  f.var <- fb$forward
  b.var <- fb$backward
  
  p.xy <- E$M
  q.x <- E$X
  q.y <- E$Y
  
  switch(j,
         "M" = e.j <- p.xy[O1[u+1], O2[v+1]],
         "X" = e.j <- q.x[1, O1[u+1]],
         "Y" = e.j <- q.y[1, O2[v+1]]
  )
  
  switch(j,
         "M" = b.j <- b.var[[j]][u+1, v+1],
         "X" = b.j <- b.var[[j]][u+1, v],
         "Y" = b.j <- b.var[[j]][u, v+1]
  )
  
  num <- f.var[[i]][u, v] * A[i, j] * e.j * b.j
  
  den <- 0
  for (k in 1:3) {
    den <- den + f.var[[k]][u, v] * b.var[[k]][u, v]
  }
  
  if (den == 0) {
    return(0)
  } else {
    xi <- num / den
    return(xi)
  }
  
}

ExeXi <- function(O, fb, A, E) {
  # Computes the new transition parameters about word w.
  #
  # Args:
  #   u, v: The values of positions on both of the observation sequences.
  #   i, j: The values of states.
  #   O: The list of the observation sequences and the lengths.
  #   fb: The list of the forward-backward algorithm.
  #   A: The matrix of transition probability.
  #   E: The list of the emission probability matrices.
  #
  # Returns:
  #   The value of the new parameters about word w.
  O1 <- O$O1
  O2 <- O$O2
  U <- O$U
  V <- O$V
  
  I <- c("M", "M", "X", "X", "X")
  J <- c("M", "X", "M", "X", "Y")
  
  for (s in 1:length(I)) {
    
    i <- I[s]
    j <- J[s]
    
    if (j=="M") {
      U_ <- U-1
      V_ <- V-1
    }
    else if (j=="X") {
      U_ <- U-1
      V_ <- V
    } else {
      U_ <- U
      V_ <- V-1
    }
    
    num <- 0
    for (u in 2:U_) {
      for (v in 2:V_) {
        num <- num + Xi(u, v, i, j, O, fb, A, E)
      }
    }
    
    den <- 0
    for (u in 2:U_) {
      for (v in 2:V_) {
        for (j_ in J) {
          den <- den + Xi(u, v, i, j_, O, fb, A, E)
        }
      }
    }
    
    a.ij_ <- num / den
    if (is.na(a.ij_)) {
      a.ij_ <- 0
    } 
    else if (is.infinite(a.ij_)) {
      a.ij_<- 0
    } else {
      # no operation
    }
    
    if (i=="M" && j=="M") {
      m.m <- a.ij_      
    }
    else if (i=="X" && j=="M") {
      xy.m <- a.ij_
    }
    else if (i=="M" && j=="X") {
      delta_ <- a.ij_
    }
    else if (i=="X" && j=="X") {
      epsilon_ <- a.ij_
    }
    else if (i=="X" && j=="Y") {
      lambda_ <- a.ij_
    }
    
  }
  
  tau.M_ <- 1-2*delta_-m.m
  tau.XY_ <- 1-epsilon_-lambda_-xy.m
  
  params1 <- c(1-2*delta_-tau.M_, 2*delta_, tau.M_)
  params1 <- params1 / sum(params1)
  params1 <- params1[-1]
  
  params2 <- c(1-epsilon_-lambda_-tau.XY_, epsilon_, lambda_, tau.XY_)
  params2 <- params2 / sum(params2)
  params2 <- params2[-1]
  
  params_ <- c(params1, params2)
  names(params_) <- c(params.name)
  params_["delta"] <- params_["delta"] / 2
  
  return(params_)
}

