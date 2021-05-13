This repository contains scripts to calculate denuation rates in landscapes where
minerals with different solubilities are present. In particular, the codes calculate
denudation rates for landscpaes where cosmogenic nuclides where measured on two
target minerals, with one being insoluble and the other one soluble.

The codes are on an adpation of Riebe and Granger equation 14. For most of the below 
scripts CronusCalc needs to be in the matlab search path.

- RiebeGranger_CarbBias: This script plots the bias in denudation rate measurements due to
enrichment/depletion of minerals for different soil depths. It utilizes exponential functions
for production rates at depth and is currently only written for 10Be.

- RiebeGranger_Cronus CarbBias: This script does the same as the above script but uses
CronusCalc to calculate production rates and is adapted for 10Be and 36Cl nuclides. 

- RiebeGranger_Cronus_PairedCRN: This is a script that computes the "real denudation rate"
given a paired nuclide measurement and known bedrock chemistry. This script computes
the parameter inversion with a neighborhood algorithm and currently does not work very well.

- RiebeGranger_Cronus_PairedCRN_MCMC:
Assuming you have a paired 10be and 36Cl measurement and know your bedrock chemistry 
(fraction of quartz in bedrock, fraciton of calcite in bedrock, fraction of other insoluble minerals
in bedrock), you can calculate the "real denudation rate". The code uses a Markov Chain Monte
Carlo approach to solve for the correct denudation rate and enrichment/depletion of minerals
in the soil.
The code also estimates the weatheirng rate based on the given parameters.

The code can be tested by running the RiebeGranger_Cronus_PairedCRN_MCMC_destData script

Richard Ott, 2021