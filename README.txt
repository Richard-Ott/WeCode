----- CDF BRANCH ---------
This repository contains functions to calculate denuation rates in landscapes where
minerals with different solubilities are present. All necessary functions are in the
subroutines folder. The repository also contains example scripts on how to use the
supplied functoins. 
In particular, the codes calculate denudation rates for landscpaes where cosmogenic 
nuclides were measured and weathering is non negligable. I have adopted the code package
for 10Be measurements on an insoluble target mineral (e.q., quartz) AND/OR 36Cl measurements
on a soluble target mineral (e.g., calcite). 
Since, this branch focuses on 10Be applications, it runs WITHOUT DECAY!!!

The codes are based on an adpation of Riebe and Granger equation 14.
Cronus Calc v2.1 is used for the calculation of production rates and needs to be in the 
Matlab search path. If you want to run the code with alluvial basin data, you also need 
TopoToolbox in your Matlab search path. he code will then compute average nuclide production
rates for your basin in a pixel-by-pixel approach. The code also provides a set of modified
Cronus Calc v2.1 functions to achieve this.

Example input files are provided for all scripts, inlcuding test data that you can use to
run the code. 

- RiebeGranger_Cronus_SingleCRN: This code calculates the "real" denudation rate from a
single nuclide measurement (soluble or insoluble target mineral), given that either a CDF
is provided. It uses a built-in optimization algorithm to find the denudation rate.

Cite as:
Ott, Richard (2022): WeCode - Weathering Corrections for denudation rates. V. 1.0. GFZ Data Services. 
https://doi.org/10.5880/GFZ.4.6.2022.001

Related references:
Marrero, S. M., Phillips, F. M., Borchers, B., Lifton, N., Aumer, R., & Balco, G. (2016). Cosmogenic 
nuclide systematics and the CRONUScalc program. Quaternary Geochronology, 31, 160–187. 
https://doi.org/10.1016/j.quageo.2015.09.005