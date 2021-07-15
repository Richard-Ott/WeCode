This repository contains scripts to calculate denuation rates in landscapes where
minerals with different solubilities are present. In particular, the codes calculate
denudation rates for landscpaes where cosmogenic nuclides where measured on one or two
target minerals, with one being insoluble and the other one soluble.

The codes are on an adpation of Riebe and Granger equation 14.
Cronus Calc v2.1 is used for the calculation of production rates and needs to be in the 
Matlab search path. If you want to run the code with alluvial basin data, you also need 
TopoToolbox in your Matlab search path.

Example input files are provided for all scripts, inlcuding test data that you can use to
run the code.

- RiebeGranger_Cronus_CarbBias: This script plots the bias in denudation rate measurements due to
enrichment/depletion of minerals for different soil depths. 

- RiebeGranger_Cronus_PairedCRN:
Assuming you have a paired 10Be and 36Cl measurement and know your bedrock or soil chemistry 
(fraction of quartz, fraction of calcite, fraction of other insoluble minerals), you 
can calculate the "real denudation rate". The code uses a Markov Chain Monte Carlo 
approach to solve for the correct denudation rate and enrichment/depletion of minerals
in the soil.
The code also estimates the weathering rate based on the given parameters.

- RiebeGranger_Cronus_SingleCRN: This code calculates the "real" denudation rate from a
single nuclide measurement (soluble or insoluble target mineral), given that either the 
bedrock or soil chemistry is provided, and the overall weathering rate is known. It also
uses a Markov Chain Monte Carlo algorithm to efficiently solve the system of equations.

Richard Ott, 2021