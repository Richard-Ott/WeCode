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

% load data
[num,sampName,X,DEMdata] = CosmoDataRead('Test_Input_Paired.xlsx',1);

%% assign data and initial basin calculations --------------------------- %

% in case your denudation rate is from an alluvial sample and you desire a
% pixel-by-pixel calculated production rate provide a DEM
% THIS STEP REQUIRES TOPOTOOLBOX FUNCTIONS (Schwanghart & Scherler, 2014)
if strcmpi('basin',DEMdata.method)
    DEMdata.DEM = GRIDobj();      % interactively choose the DEM that encompasses the basin
    DEMdata.export = 1;           % do you want to save the individual sample scaling factors as .mat file?
    % This can be useful when the computation for scaling schemes like 'sa'
    % and 'sf'  takes a long time for a big basin and you want the scaling
    % factors saved for later
    
    [DEMdata.DB,DEMdata.utmzone] = getBasins(DEMdata.DEM,num.num10(:,2),num.num10(:,1),'ll');  % delineate drainage basins and check their geometry
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

% Report the values %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['Denudation rate = ' num2str(round(MAP)) ' ' char(177) ' ' num2str(round(MAP_uncerts)) ' mm/ka'])
disp(['The calculated weathering rate = ' num2str(round(WMAP)) ' ' char(177) ' ' num2str(round(W_uncert)) ' mm/ka'])
switch X.mode
    case 'soil'
        disp(['Fraction of quartz in bedrock fQzB = '  num2str(round(XMAP.fQzB,2)) ' ' char(177) ' ' num2str(round(X_uncerts(1),2))])
        disp(['Fraction of X in bedrock fXB = '        num2str(round(XMAP.fXB,2))  ' ' char(177) ' ' num2str(round(X_uncerts(2),2))])  
        disp(['Fraction of calcite in bedrock fCaB = ' num2str(round(XMAP.fCaB,2)) ' ' char(177) ' ' num2str(round(X_uncerts(3),2))])     
    case 'bedrock'
        disp(['Fraction of quartz in soil fQzS = '     num2str(round(XMAP.fQzS,2)) ' ' char(177) ' ' num2str(round(X_uncerts(1),2)) ])
        disp(['Fraction of X in soil fXS = '           num2str(round(XMAP.fXS,2))  ' ' char(177) ' ' num2str(round(X_uncerts(2),2))])
        disp(['Fraction of calcite in soil fCaS = '    num2str(round(XMAP.fCaS,2)) ' ' char(177) ' ' num2str(round(X_uncerts(3),2))])
end

% export = input('Do you want to export your results? "y" or "n"? ','s');
% if strcmpi(export,'y')
%     vars = {'Name','fQzS','fCaS','fXS','fQzB','fCaB','fXB','W','D'};
%     out_table = table(txt10,X.fQzS,X.fCaS,X.fXS,X.fQzB,X.fCaB,X.fXB,W,MAP(2),'VariableNames',vars);
%     writetable(out_table,[ '.\output\JO\pairedN\' txt10{1} '_' X.mode '.xlsx'])
% end
