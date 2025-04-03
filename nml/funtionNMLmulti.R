# install.packages("devtools")
# run the next line if you already have rstan installed
# remove.packages(c("StanHeaders", "rstan"))
# install.packages("devtools")
# devtools::install_github("davidkellen/NMLmulti")
remotes::install_github("davidkellen/NMLmulti")
library(NMLmulti)
nml_scr <- run_nml(fun= SCR, parl=10, ks=rep(3,5), Ns=c(rep(100,4),200), 
                   fits = 2, batchsize=5000, burn=10000, precision=0.1)

# ___________________________________________________________________________ #

Rcpp::Rcpp.package.skeleton("examplemodels", attributes = TRUE, 
                            cpp_files = "examplemodels.cpp")
devtools::install("examplemodels")

library("NMLmulti")
library("examplemodels")
remove.packages("examplemodels")
# remove.packages("NMLmulti")
nml_scr <- run_nml(fun= SCR2, parl=10, ks=rep(3,5), Ns=c(rep(100,4),200),
                   packages_multicore = "examplemodels",
                   fits = 2, batchsize=5000, burn=10000, precision=0.1)

library(NMLmulti)
## basic example code

