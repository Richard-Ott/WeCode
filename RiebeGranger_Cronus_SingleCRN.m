% The script shows how to calculate the denudation rate from a single nuclide
% measurement of an insoluble target mineral with a CDF value.
%
% Richard Ott, 2021
%
% the current version is written for 10Be and 36Cl, could easily be
% expanded to other nuclides

clc
clear 
close all
addpath '.\subroutines'

% load data
[num,sampName,X,DEMdata] = CosmoDataRead('Test_Input_Single_Be.xlsx');

%% assign data and initial basin calculations --------------------------- %

% in case your denudation rate is from an alluvial sample and you desire a
% pixel-by-pixel calculated production rate provide a DEM
% THIS STEP REQUIRES TOPOTOOLBOX FUNCTIONS (Schwanghart & Scherler, 2014)
if strcmpi('basin',DEMdata.method)
    DEMdata.DEM = GRIDobj();        % interactively choose the DEM that encompasses the basin
    DEMdata.export = 1;             % do you want to save the individual sample scaling factors as .mat file?
    % This can be useful when the computation for scaling schemes like 'sa'
    % and 'sf'  takes a long time for a big basin and you want the scaling
    % factors saved for later
    
    [DEMdata.DB,DEMdata.utmzone] = getBasins(DEMdata.DEM,num(:,2),num(:,1),'ll');  % delineate drainage basins and check their geometry
end


%% Calculate production rates ------------------------------------------- %

Cronus_prep = {@Cronus_prep10, @Cronus_prep36};

pars = Cronus_prep{X.n}(num,DEMdata);

%% Run Optimization for "real" denudation rate -------------------------- %

D = [50,1e3];                          % Denudation min/max in mm/ka (Dmin > Weathering rate)
thres = 0.1;                           % threshold of nuclide concentration error for optimization in % of N

% run inversion
MAP = singleCRN_Optim(pars,D,X,thres);   

% OPTIONAL
% estimate uncertainty, this take quite some time to calculate, therefore I
% commented the next line to speed up the example
MAP_uncerts = singleCRN_uncerts(pars,D,X,MAP,thres);

%% OUTPUT RESULTS ------------------------------------------------------- %

% outputs the values for the higher rate if there was a nonunique solution
try
    % disp(['Denudation rate = ' num2str(round(MAP(end))) ' mm/ka'])
    disp(['Denudation rate = ' num2str(round(MAP(end))) ' ' char(177) ' ' num2str(round(MAP_uncerts)) ' mm/ka'])
end

