This repository contains scripts to calculate denuation rates in landscapes where
minerals with different solubilities are present. In particular, the codes calculate
denudation rates for landscpaes where cosmogenic nuclides where measured on one or two
target minerals, with one being insoluble and the other one soluble.

The codes are on an adpation of Riebe and Granger equation 14. For most of the below 
scripts CronusCalc needs to be in the matlab search path.

- RiebeGranger_CarbBias: This script plots the bias in denudation rate measurements due to
enrichment/depletion of minerals for different soil depths. It utilizes exponential functions
for production rates at depth and is currently only written for 10Be.

- RiebeGranger_Cronus CarbBias: This script does the same as the above script but uses
CronusCalc to calculate production rates and is adapted for 10Be and 36Cl nuclides. 

- RiebeGranger_Cronus_PairedCRN_NA (NOT RECOMMENDED): This is a script that computes the "real 
enudation rate" given a paired nuclide measurement and known bedrock or soil chemistry. This 
script computes the parameter inversion with a neighborhood algorithm and IS OUTDATED AND DOES 
NOT WORK WELL.

- RiebeGranger_Cronus_PairedCRN:
Assuming you have a paired 10Be and 36Cl measurement and know your bedrock or soil chemistry 
(fraction of quartz, fraction of calcite, fraction of other insoluble minerals), you 
can calculate the "real denudation rate". The code uses a Markov Chain Monte Carlo 
approach to solve for the correct denudation rate and enrichment/depletion of minerals
in the soil.
The code also estimates the weatheirng rate based on the given parameters.
You can test this code with the test data set provided. Just set test = 1 in the script.

- RiebeGranger_Cronus_SingleCRN: This code calculates the "real" denudation rate from a
single nuclide measurement (soluble or insoluble target mineral), given that either the 
bedrock or soil chemistry is provided, and the overall weathering rate is known. It also
uses a Markov Chain Monte Carlo algorithm to efficiently solve the system of equations.
You can test this code with the test data set provided.

Richard Ott, 2021