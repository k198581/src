library(doMC)

#cores <- detectCores() / 2
cores <- detectCores()
registerDoMC(cores)
