% This in an example script on how to calculate the denudation rate for a 
% paired nuclide measurement of a soluble and an insoluble target mineral. 
%
% The parameter search is performed via a Markov-Chain Monte Carlo approach
% with a Metropolis Hastings sampling alrogithm. 
%
% Richard Ott, 2021

% the current version is written for 10Be and 36Cl, could easily be
% expanded to other nuclides
clc
clear 
close all
addpath '.\subroutines'

Dout = nan(1,5); D_uncerts = nan(size(Dout)); Wout = nan(size(Dout)); Wout_uncerts = nan(size(Dout));
Xout= cell(1,5);
Xout_uncert = cell(1,5);
for i = 2:5
% load data
[num,sampName,X,DEMdata] = CosmoDataRead('pairedCRN_Cronus_Basin.xlsx',i);

%% assign data and initial basin calculations --------------------------- %

% in case your denudation rate is from an alluvial sample and you desire a
% pixel-by-pixel calculated production rate provide a DEM
% THIS STEP REQUIRES TOPOTOOLBOX FUNCTIONS (Schwanghart & Scherler, 2014)
if strcmpi('basin',DEMdata.method)
    DEMdata.DEM = GRIDobj('C:\Users\r_ott\Dropbox\Richard\NAGRA\GIS\SRTM30\SRTM30_DEM_utm.tif');      % interactively choose the DEM that encompasses the basin
    DEMdata.export = 1;           % do you want to save the individual sample scaling factors as .mat file?
    % This can be useful when the computation for scaling schemes like 'sa'
    % and 'sf'  takes a long time for a big basin and you want the scaling
    % factors saved for later
    
    [DEMdata.DB,DEMdata.utmzone] = getBasins_nonInteractive(DEMdata.DEM,num.num10(:,2),num.num10(:,1),'ll',32);  % delineate drainage basins and check their geometry
end


%% Calculate production rates ------------------------------------------- %

pars10 = Cronus_prep10(num.num10,DEMdata);
if isnan(num.num36(11)); num.num36(11) = pars10.nominal10(13); end % if the 10Be effective attenuation length was already computed, don't do the same again for 36Cl (saves time)
pars36 = Cronus_prep36(num.num36,DEMdata);


%% Run MCMC inversion for "real" denudation rate and weathering rate ---- %

D = [50,5e2];       % Denudation rates to test, min/max in mm/ka 

% run inversion
[XMAP,MAP,WMAP] = pairedCRN_Optim(pars10,pars36,D,X);

% OPTIONAL %%%%%%%%
% estimate uncertainty, this take quite some time to calculate, therefore I
% commented the next line to speed up the example
[MAP_uncerts, X_uncerts, W_uncert] = pairedCRN_uncerts(pars10,pars36,D,X,XMAP,MAP,WMAP);

%% REPORT RESULTS ------------------------------------------------------- %

% store the values %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Dout(i)  = MAP;
D_uncerts(i) = MAP_uncerts;
Wout(i) = WMAP;
Wout_uncerts(i) = W_uncert;
Xout{i} = XMAP;
Xout_uncert{i} = X_uncerts;

if i == 1
    vars = {'Name','fQzS','dfQzS','fXS','dfXS','fCaS','dfCaS','fQzB','fCaB','fXB','W','dW','D','dD'};
    out_table = table({sampName},XMAP.fQzS,X_uncerts(1),XMAP.fXS,X_uncerts(2),XMAP.fCaS,X_uncerts(3),...
        XMAP.fQzB,XMAP.fCaB,XMAP.fXB, Wout(i),Wout_uncerts(i), Dout(i), D_uncerts(i) ,'VariableNames',vars);
else
    out_table = [out_table; table({sampName},XMAP.fQzS,X_uncerts(1),XMAP.fXS,X_uncerts(2),XMAP.fCaS,X_uncerts(3),...
        XMAP.fQzB,XMAP.fCaB,XMAP.fXB, Wout(i),Wout_uncerts(i), Dout(i), D_uncerts(i) ,'VariableNames',vars)];
end

end
writetable(out_table,[ '.\output\JO\pairedN\PairedNuclide_' X.mode '.xlsx'])

