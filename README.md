# A-J Multinomial Process Tree Model

## Overview

This repository contains scripts and data for evaluating recognition models using the Extended-MPT framework, based on the original experimental data from Juola et al. (2019). Data: https://osf.io/y78mk/.

## Repository Contents

Scripts (place these and the main data file in the same directory)
- aj.mpt.R — Main script for model fitting, data analysis, and figure generation. Automatically calls required secondary scripts.
- aj.sim.R — Executes simulations of multinomial models. Automatically calls required secondary scripts.
- aj.fun.R — Auxiliary functions used by aj.mpt.R and aj.sim.R.

Folder Structure (first level)
- Simulation Performance Measures
  * pm.R — simulations and performance measures (results of these simulations).
  * nml2ht.RData, nmlsdt.RData, nmlaj.RData — files underlying Appendix C; produced by running aj.sim.R.
- model_files/ — TXT model specifications.
- figures/ — figures and tables created by aj.mpt.R.
- main_fits/ — fit files created by aj.mpt.R.
- nml/ — implementation and computation of normalized maximum likelihood (NML) penalties (C++). Relevant sources: model2ht.cpp, modelaj.cpp, modelsdt.cpp. For NML details, see Kellen & Klauer (2020). These files are required to generate nml2ht.RData, nmlsdt.RData, and nmlaj.RData.
- data/ — recategorized data used in the analyses: d.data.cl.rt.RData. Supplementary recategorized data for five-condition manipulations: data.5j.2CL.RData and data.5j.3CL.RData. For the original, non-transformed data see https://osf.io/y78mk/ (file d.juola).

## Usage

1) Import the repository: use the entire folder structure.
2) Execute results: run aj.mpt.R to produce main fit files (saved in main_fits) and visualization files (saved in figures).
3) Simulations results: pm.R contains the performance results for the simulations in Appendix C. To replicate these simulations, run aj.sim.R. Both pm.R and the NML penalty files (nml2ht.RData, nmlsdt.RData, nmlaj.RData) are created by executing aj.sim.R.
4) Functions, models, and penalties: ensure aj.fun.R, the model_files folder, and the nml folder are in the same directory as aj.sim.R and aj.mpt.R; they reference each other.


## Reference

Juola, J. F., Caballero-Sanz, A., Muñoz-García, A. R., Botella, J., & Suero, M. (2019). Familiarity, recollection, and receiver-operating characteristic (ROC) curves in recognition memory. Memory & Cognition, 47(4), 855–876.
Kellen, D., & Klauer, K. C. (2020). Selecting amongst multinomial models: An apologia for normalized maximum likelihood. Journal of Mathematical Psychology, 97, 1–47. https://doi.org/10.1016/j.jmp.2020.102367
