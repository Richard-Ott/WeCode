% This script shows how to calculate denudation rates in the case of
% soil-bedrock weathering.
% Richard Ott, 2021

clc
clear 
close all
addpath '.\subroutines'

% load data
[num,sampName,X,DEMdata] = CosmoDataRead('Test_Input_Single_Be.xlsx');

%% assign data and initial basin calculations --------------------------- %

% in case your denudation rate is from an alluvial sample and you desire a
% pixel-by-pixel calculated production rate, please provide a DEM
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

%% Calculate corrected denudation rate

D = SBW_search(pars,X);  % denudation rate in mm/ka

% calculate uncertinty from N,W, and productioon rates
D_uncerts = SBW_uncerts(pars,X,D);

disp(['The corrected denudation rate is ' num2str(round(D)) ' ' char(177) ' ' num2str(round(D_uncerts))  ' mm/ka'])
