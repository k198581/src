source("phmm/initialization.R")
source("phmm/forward.R")
source("phmm/backward.R")

Po.vec <- c(NULL)
for (i in 1:100) {
  print(i)
  source("phmm/initialization.R")
  a <- Forward(O, params, E)
  b <- Backward(O, params, E)
  Po.vec <- c(Po.vec, abs(a - b))
}
