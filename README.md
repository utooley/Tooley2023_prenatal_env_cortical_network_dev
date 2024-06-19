# Tooley2023_prenatal_env_cortical_network_dev

**Author**: Ursula Tooley, 2023-2024

This repository contains code for the analyses described in our paper:

"Prenatal environment is associated with the pace of cortical network development over the first three years of life.""

Preprint available [here](https://doi.org/10.1101/2023.08.18.552639).

Final paper link will be available here.

# Dependencies & operating system
This code was run on Mac OS Monterey 12.6.3.

*MATLAB version 2021b* and *R* version 4.3.1 (2023-06-16) were used to conduct analyses in this paper. 

Preprocessing was performed with [4dfp](https://4dfp.readthedocs.io/en/latest/), Melbourne Children’s Regional Brain Atlases (MCRIB), Connectome Workbench 1.2.3, and Freesurfer 7.2.

# Organization

The *gamm_models* folder contains a script with functions necessary for running the code called *gamm_functions.R*. This script is sourced in the main analysis file, titled `manuscript_analyses.R` and the supplementary analysis file, titled `sensitivity_analyses.R`, prior to continuing analyses.
