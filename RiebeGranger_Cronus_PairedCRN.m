% The script calculates the denudation rate for a paired nuclide
% measurement of a soluble and an insoluble target mineral. 
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

% Priors
D = [5,1e3];                           % Denudation min/max in mm/ka
% second prior would be either bedrock or soil quartz content and is
% automatically set to fractions between 0 and 1 in the MCMC function.

% run inversion
[X,MAP,post,W] = pairedCRN_MCMC(pars10,pars36,D,X);


%% OUTPUT RESULTS ------------------------------------------------------- %

% plot the MCMC chains 
figure()
subplot(1,3,1); plot(post(:,1))
xlabel('iteration'); ylabel('fraction Qz'); ylim([0,1])
subplot(1,3,2); plot(post(:,2))
xlabel('iteration'); ylabel('Denudation rate g/cm²/ka'); ylim(D)
subplot(1,3,3); plot(post(:,3))
xlabel('iteration'); ylabel('log posterior probability')

% Report the values
disp(['Denudation rate = ' num2str(round(MAP(2))) ' mm/ka'])
switch X.mode
    case 'soil'
        disp(['Fraction of quartz in bedrock fQzB = ' num2str(X.fQzB)])
        disp(['Fraction of X in bedrock fXB = ' num2str(X.fXB)])  
        disp(['Fraction of calcite in bedrock fCaB = ' num2str(X.fCaB)])        
    case 'bedrock'
        disp(['Fraction of quartz in soil fQzS = ' num2str(MAP(1))])
        disp(['Fraction of X in soil fXS = ' num2str(X.fXS)])
        disp(['Fraction of calcite in soil fCaS = ' num2str(X.fCaS)])
end
disp(['The calculated weathering rate = ' num2str(round(W)) ' mm/ka'])

export = input('Do you want to export your results? "y" or "n"? ','s');
if strcmpi(export,'y')
    vars = {'Name','fQzS','fCaS','fXS','fQzB','fCaB','fXB','W','D'};
    out_table = table(txt10,X.fQzS,X.fCaS,X.fXS,X.fQzB,X.fCaB,X.fXB,W,MAP(2),'VariableNames',vars);
    writetable(out_table,[ '.\output\JO\pairedN\' txt10{1} '_' X.mode '.xlsx'])
end
