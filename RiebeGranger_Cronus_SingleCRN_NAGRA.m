% The script shows how to calculate the denudation rate from a single nuclide
% measurement of a soluble or an insoluble target mineral. 
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

Dout = nan(1,5); D_uncerts = nan(size(Dout)); 
Xout= cell(1,5);
Xout_uncert = cell(1,5);
utmzone = 32;

for i= 1:10
% load data
[num,sampName,X,DEMdata] = CosmoDataRead('10Be_Cronus_Basin_Soil.xlsx',i);

%% assign data and initial basin calculations --------------------------- %

% in case your denudation rate is from an alluvial sample and you desire a
% pixel-by-pixel calculated production rate provide a DEM
% THIS STEP REQUIRES TOPOTOOLBOX FUNCTIONS (Schwanghart & Scherler, 2014)
if strcmpi('basin',DEMdata.method)
    DEMdata.DEM = GRIDobj('C:\Users\r_ott\Dropbox\Richard\NAGRA\GIS\SRTM30\SRTM30_DEM_utm.tif');        % interactively choose the DEM that encompasses the basin
    DEMdata.export = 1;             % do you want to save the individual sample scaling factors as .mat file?
    % This can be useful when the computation for scaling schemes like 'sa'
    % and 'sf'  takes a long time for a big basin and you want the scaling
    % factors saved for later
    
    [DEMdata.DB,DEMdata.utmzone] = getBasins_nonInteractive(DEMdata.DEM,num(:,2),num(:,1),'ll',utmzone);  % delineate drainage basins and check their geometry
end


%% Calculate production rates ------------------------------------------- %

Cronus_prep = {@Cronus_prep10, @Cronus_prep36};

pars = Cronus_prep{X.n}(num,DEMdata);

%% Run Optimization for "real" denudation rate -------------------------- %
% PRIORS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
D = [50,5e2];                          % Denudation min/max in mm/ka (Dmin > Weathering rate)
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

% OPTIONAL
% estimate uncertainty, this take quite some time to calculate, therefore I
% commented the next line to speed up the example
[MAP_uncerts, X_uncerts] = singleCRN_uncerts(pars,D,X,X_MAP,MAP,thres);

%% OUTPUT RESULTS ------------------------------------------------------- %

% outputs the values for the higher rate if there was a nonunique solution
% store the values %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Dout(i)  = MAP;
D_uncerts(i) = MAP_uncerts;
Xout{i} = X_MAP;
Xout_uncert{i} = X_uncerts;

if i == 1
    vars = {'Name','fQzS','dfQzS','fXS','dfXS','fCaS','dfCaS','fQzB','fXB','fCaB','D','dD'};
    out_table = table({sampName},X_MAP.fQzS,X_uncerts(1),X_MAP.fXS,X_uncerts(2),X_MAP.fCaS,X_uncerts(3),...
        X_MAP.fQzB,X_MAP.fCaB,X_MAP.fXB, Dout(i), D_uncerts(i) ,'VariableNames',vars);
else
    out_table = [out_table; table({sampName},X_MAP.fQzS,X_uncerts(1),X_MAP.fXS,X_uncerts(2),X_MAP.fCaS,X_uncerts(3),...
        X_MAP.fQzB,X_MAP.fCaB,X_MAP.fXB, Dout(i), D_uncerts(i) ,'VariableNames',vars)];
end
disp(i)
end

writetable(out_table,[ '.\output\JO\10Be\Be10_' X.mode '.xlsx'])
