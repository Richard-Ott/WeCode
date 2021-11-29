This repository contains functions to calculate denuation rates in landscapes where
minerals with different solubilities are present. All necessary functions are in the
subroutines folder. The repository also contains example scripts on how to use the
supplied functoins. 
In particular, the codes calculate denudation rates for landscpaes where cosmogenic 
nuclides were measured and weathering is non negligable. I have adopted the code package
for 10Be measurements on an insoluble target mineral (e.q., quartz) AND/OR 36Cl measurements
on a soluble target mineral (e.g., calcite). 

The codes are based on an adpation of Riebe and Granger equation 14.
Cronus Calc v2.1 is used for the calculation of production rates and needs to be in the 
Matlab search path. If you want to run the code with alluvial basin data, you also need 
TopoToolbox in your Matlab search path. The code will then compute average nuclide production
rates for your basin in a pixel-by-pixel approach. The code also provides a set of modified
Cronus Calc v2.1 functions. 

Example input files are provided for all scripts, inlcuding test data that you can use to
run the code. 
Single nuclide measurement with an independent weathering rate:
Test_Input_Single.xlsx - >10Be, Test_Input_Single2.xlsx -> 36Cl
Test_input_Paired.xlsx is for running the paired nuclide inversion

- RiebeGranger_Cronus_PairedCRN:
Assuming you have a paired 10Be and 36Cl measurement and know your bedrock OR soil chemistry 
(fraction of quartz, fraction of calcite, fraction of other insoluble minerals), you 
can calculate the "real denudation rate". The code uses an optimization algorithm to solve 
for the correct denudation rate and enrichment/depletion of minerals in the soil.
The code also estimates the weathering rate based on the given parameters.

- RiebeGranger_Cronus_SingleCRN: This code calculates the "real" denudation rate from a
single nuclide measurement (soluble or insoluble target mineral), given that either the 
bedrock or soil chemistry is provided, and the overall weathering rate is known. It uses
a built-in optimization algorithm to find the denudation rate.

- Soil_Bedrock_weathering: This script calculates denudation rates for cases where weathering
concentrates along the soil-bedrock interface. To solve the equations it uses a the fminsearch
Matlab optimization algorithm. Run this script with one of the single nuclide test files.

- RiebeGranger_Cronus_CarbBias: This script plots the bias in denudation rate measurements due to
enrichment/depletion of minerals for different soil depths. 

Richard Ott, 2021