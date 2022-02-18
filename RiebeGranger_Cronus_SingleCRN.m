% The script shows how to calculate the denudation rate from a single nuclide
% measurement of a soluble or an insoluble target mineral. Use the provided
% single CRN input sheets.
%
% The parameter search is performed via a Markov-Chain Monte Carlo approach
% with a Metropolis Hastings sampling alrogithm. 
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
[num,sampName,X,DEMdata] = CosmoDataRead('Test_Input_Single_Cl.xlsx');

%% assign data and initial basin calculations --------------------------- %

% in case your denudation rate is from an alluvial sample and you desire a
% pixel-by-pixel calculated production rate provide a DEM
% THIS STEP REQUIRES TOPOTOOLBOX FUNCTIONS (Schwanghart & Scherler, 2014)
if strcmpi('basin',DEMdata.method)
    DEMdata.DEM = GRIDobj();        % interactively choose the DEM that encompasses the basin
    % Scaling schemes like 'sa' and 'sf'  can take a long time for a big 
    % basin and you want the scaling. You may want to save the scaling data
    % for later re-runs.
    
    [DEMdata.DB,DEMdata.utmzone] = getBasins(DEMdata.DEM,num(:,2),num(:,1),'ll');  % delineate drainage basins and check their geometry
end


%% Calculate production rates ------------------------------------------- %

Cronus_prep = {@Cronus_prep10, @Cronus_prep36};

pars = Cronus_prep{X.n}(num,DEMdata);

%% Run Optimization for "real" denudation rate -------------------------- %
% PRIORS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
D = [50,1e3];                          % Denudation min/max in mm/ka (Dmin > Weathering rate)
thres = 0.1;                           % threshold of nuclide concentration error at which optimization has converged in % of N

% if you're looking at the soluble target mineral, there might be more than
% one solution. This function will check for you if the parameters you
% suppied might lead to a nonunique solution. It's not necessary to run
% this, but advisable.
if X.n == 2 
    pars = solCRN_check(pars,D,X,thres,0);
end

% run inversion
[MAP,X_MAP] = singleCRN_Optim(pars,D,X,thres);   
% MAP = singleCRN_Optim(pars,D,X,thres); % run this in case of 10Be and no composition


% OPTIONAL
% estimate uncertainty, this take quite some time to calculate, therefore I
% commented the next line to speed up the example
[MAP_uncerts, X_uncerts] = singleCRN_uncerts(pars,D,X,MAP,thres,X_MAP);
% MAP_uncerts = singleCRN_uncerts(pars,D,X,MAP,thres); % run this in case of 10Be and no composition
%% OUTPUT RESULTS ------------------------------------------------------- %

% outputs the values for the higher rate if there was a nonunique solution
try
    disp(['Denudation rate = ' num2str(round(MAP(end))) ' ' char(177) ' ' num2str(round(MAP_uncerts)) ' mm/ka'])
    switch X.mode
        case 'soil'
            disp(['Fraction of quartz in bedrock fQzB = '  num2str(round(X_MAP(end).fQzB,2)) ' ' char(177) ' ' num2str(round(X_uncerts(1),2))])
            disp(['Fraction of X in bedrock fXB = '        num2str(round(X_MAP(end).fXB,2))  ' ' char(177) ' ' num2str(round(X_uncerts(2),2))])  
            disp(['Fraction of calcite in bedrock fCaB = ' num2str(round(X_MAP(end).fCaB,2)) ' ' char(177) ' ' num2str(round(X_uncerts(3),2))])     
        case 'bedrock'
            disp(['Fraction of quartz in soil fQzS = '     num2str(round(X_MAP(end).fQzS,2)) ' ' char(177) ' ' num2str(round(X_uncerts(1),2))])
            disp(['Fraction of X in soil fXS = '           num2str(round(X_MAP(end).fXS,2))  ' ' char(177) ' ' num2str(round(X_uncerts(2),2))])
            disp(['Fraction of calcite in soil fCaS = '    num2str(round(X_MAP(end).fCaS,2)) ' ' char(177) ' ' num2str(round(X_uncerts(3),2))])
    end
end

