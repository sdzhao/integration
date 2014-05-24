## ##############################################################
## README, 12.3.2013
## More powerful genetic association testing via a new statistical framework for integrative genomics
## Sihai D. Zhao, T. Tony Cai, and Hongzhe Li
## ##############################################################
## ==============================================================
## Description
## ==============================================================
R implementation of methods and simulations. Each simulation file was designed to run on a a computing cluster. Each file contains an "env" variable, which serves as the random seed for the data generating process and can be set as an environment variable using the cluster job submission system.

## ==============================================================
## Contents
## ==============================================================
direct.R
- direct effect; Example 2 in main text

easy.R
- easy setting, model is correctly specified; Example 1 in main text

fxns.R
- implementation of proposed method for cohort and case-control sampling
- descriptions of function arguments are given before the function definition

hd.R
- high dimensional expression data; Example 6 in main text

int.R
- interactions; Example 5 of main text

measerr.R
- measurement error; Example 3 of main text

mis.R
- misspecified gene set; Example 4 of main text
