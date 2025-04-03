
# A-J Multinomial Process Tree Model

## Overview

This repository contains scripts and data for evaluating recognition models using the Extended-MPT framework, based on the original experimental data from Juola et al. (2019). The data can be accessed at: https://osf.io/y78mk/.

## Repository Contents

### Scripts
- Main Scripts: There are three primary scripts necessary to replicate our main results. These scripts and the main data file must be located in the same directory to function properly.

`aj.mpt.R`: Main script for model fitting, data analysis, and figure generation. This script automatically calls the required secondary scripts.
`aj.sim.R`: Script dedicated to running simulations of multinomial models.

Note: The main scripts and data files must be located in the same directory for proper execution.

### Data Files
- Main Data: The principal recategorized data file used in the analyses is `d.data.cl.rt.RData`.
- Supplementary Data: Recategorized data files are for analyses involving five conditions of manipulation of relative target frequencies, `data.5j.2CL.RData` and `data.5j.3CL.RData`. 

### Models and Functions
- Model Specifications: Detailed in txt files in the `model_files` folder. Some models are also implemented in C++ to optimize the computation of the normalized maximum likelihood (NML) penalty. The relevant C++ files are: `model2ht.cpp`, `modelaj.cpp` and `modelsdt.cpp`. For further explanation of NML penalties, see Kellen & Klauer,2020.
- Custom Functions: Developed specifically for this study and located in `aj.fun.R`.
- Analysis and Visualization: All operations for fitting models, analyzing data, and generating figures are performed by running `aj.mpt.R`. This script will automatically call the necessary secondary scripts.

### Output
- Fit Objects: Main fits objects are stored in the `main_fit` folder. Running the main script, `aj.mpt.R` or `aj.sim.R`, generates all required objects, including those in the `tables` and `main_fit` folders.

## Reference

Juola, J. F., Caballero-Sanz, A., Muñoz-García, A. R., Botella, J., & Suero, M. (2019). Familiarity, recollection, and receiver-operating characteristic (ROC) curves in recognition memory. Memory & Cognition, 47(4), 855-876.
Kellen, D., & Klauer, K. C. (2020). Selecting amongst multinomial models: An apologia for normalized maximum likelihood. Journal of Mathematical Psychology, 97, 1–47. https://doi.org/10.1016/j.jmp.2020.102367

