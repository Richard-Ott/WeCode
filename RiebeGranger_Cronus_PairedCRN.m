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
[num,sampName,X,DEMdata] = CosmoDataRead('Test_Input_Paired.xlsx');

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
pars36 = Cronus_prep36(num.num36,DEMdata);


%% Run MCMC inversion for "real" denudation rate and weathering rate ---- %

D = [5,1e3];       % Denudation rates to test, min/max in mm/ka

% run inversion
[XMAP,MAP,post,WMAP] = pairedCRN_MCMC(pars10,pars36,D,X);

% OPTIONAL %%%%%%%%
% estimate uncertainty, this take quite some time to calculate, therefore I
% commented the next line to speed up the example
[MAP_uncerts, X_uncerts, W_uncert] = pairedCRN_uncerts(pars10,pars36,D,X,XMAP,MAP,WMAP);

%% REPORT RESULTS ------------------------------------------------------- %

% investigate the MCMC chains (optional, just to check the MCMC convergence)
figure()
subplot(1,3,1); plot(post(:,1),'LineWidth',1.5)
xlabel('iteration'); ylabel('fraction Qz'); ylim([0,1])
subplot(1,3,2); plot(post(:,2),'LineWidth',1.5)
xlabel('iteration'); ylabel('Denudation rate mm/ka'); ylim(D)
subplot(1,3,3); plot(post(:,3),'LineWidth',1.5)
xlabel('iteration'); ylabel('log posterior probability')

% Report the values %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% disp(['Denudation rate = ' num2str(round(MAP(2))) ' mm/ka'])
% disp(['The calculated weathering rate = ' num2str(round(W)) ' mm/ka'])
disp(['Denudation rate = ' num2str(round(MAP(2))) ' ' char(177) ' ' num2str(round(MAP_uncerts)) ' mm/ka'])
disp(['The calculated weathering rate = ' num2str(round(WMAP)) ' ' char(177) ' ' num2str(round(W_uncert)) ' mm/ka'])
switch X.mode
    case 'soil'
        disp(['Fraction of quartz in bedrock fQzB = ' num2str(XMAP.fQzB) ' ' char(177) ' ' num2str(round(X_uncerts(1),3))])
        disp(['Fraction of X in bedrock fXB = ' num2str(XMAP.fXB) ' ' char(177) ' ' num2str(round(X_uncerts(2),3))])  
        disp(['Fraction of calcite in bedrock fCaB = ' num2str(XMAP.fCaB) ' ' char(177) ' ' num2str(round(X_uncerts(3),3))])     
%         disp(['Fraction of quartz in bedrock fQzB = ' num2str(X_MAP.fQzB) ])
%         disp(['Fraction of X in bedrock fXB = ' num2str(X_MAP.fXB)])  
%         disp(['Fraction of calcite in bedrock fCaB = ' num2str(X_MAP.fCaB)])  
    case 'bedrock'
        disp(['Fraction of quartz in soil fQzS = ' num2str(XMAP.fQzS) ' ' char(177) ' ' num2str(round(X_uncerts(1),3)) ])
        disp(['Fraction of X in soil fXS = ' num2str(XMAP.fXS) ' ' char(177) ' ' num2str(round(X_uncerts(2),3))])
        disp(['Fraction of calcite in soil fCaS = ' num2str(XMAP.fCaS) ' ' char(177) ' ' num2str(round(X_uncerts(3),3))])
%         disp(['Fraction of quartz in soil fQzS = ' num2str(X_MAP.fQzS) ])
%         disp(['Fraction of X in soil fXS = ' num2str(X_MAP.fXS)])
%         disp(['Fraction of calcite in soil fCaS = ' num2str(X_MAP.fCaS)])
end

% export = input('Do you want to export your results? "y" or "n"? ','s');
% if strcmpi(export,'y')
%     vars = {'Name','fQzS','fCaS','fXS','fQzB','fCaB','fXB','W','D'};
%     out_table = table(txt10,X.fQzS,X.fCaS,X.fXS,X.fQzB,X.fCaB,X.fXB,W,MAP(2),'VariableNames',vars);
%     writetable(out_table,[ '.\output\JO\pairedN\' txt10{1} '_' X.mode '.xlsx'])
% end
