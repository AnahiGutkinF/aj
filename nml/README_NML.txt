NMLmulti Model Estimation Example
=================================

This guide explains how to estimate models using the NMLmulti R package, including integrating custom models written in C++.

Requirements:
-------------
- R (version 4.0 or higher)
- R packages: NMLmulti, Rcpp

Steps:
------

1) Install and load the NMLmulti package:

   remotes::install_github("davidkellen/NMLmulti")
   library(NMLmulti)

2) Create an R package from a C++ model file (e.g., for 2HT, SDT, or AJ models):

   Rcpp::Rcpp.package.skeleton("examplemodels", 
                                attributes = TRUE, 
                                cpp_files = "examplemodels.cpp")

   Note: Make sure the .cpp file (e.g., "examplemodels.cpp") contains the model function. For replication of our study use "modelsdt.cpp", "model2ht.cpp" and "modelaj.cpp". This should create new folder named after the cpp file. For instance "modelsdt".

3) Run the NML estimation using the run_nml() function:

   nml_scr <- run_nml(
     fun = SCR2,
     parl = 10,
     ks = rep(3, 5),
     Ns = c(rep(100, 4), 200),
     packages_multicore = "examplemodels",
     fits = 2,
     batchsize = 5000,
     burn = 10000,
     precision = 0.1
   )

Additional Notes:
-----------------
- Replace SCR2 with the actual function name defined in your .cpp model.
- Make sure the package "examplemodels" is installed and available before running run_nml().
- The estimation process uses multicore processing. Adjust 'batchsize', 'burn', and 'precision' based on your data and computational resources.