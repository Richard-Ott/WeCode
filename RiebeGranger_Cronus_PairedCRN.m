% This script calculates the denudation rate from bedrock containing a mix
% of insoluble an soluble minerals.
% Loosely based on Riebe and Granger eqn. 14.
% Conversly to Riebe and Granger, 2014, the nuclide concentrations are not
% calculated with exponentials. Nuclide concentraions are calculated used
% CronusCalc (Marrero, 2016). 
% The script calculates the denudation rate for a paired nuclide
% measurement of a soluble and an insoluble target mineral. Despite, the
% nuclide cooncentrations, the bedrock OR soil chemistry is needed. This code
% assumes that bedrock and soil only consist of the two target minerals and
% another mineral that is assumed to be insoluble.
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

test = 0; % Do you want to run this inversion with the test data set to explore functionality? (0-no, 1-yes)
if test
    load test_data_input_v2.mat
else      % otherwise please load your respective data files

% addpath 'C:\Users\r_ott\Dropbox\Richard\Crete\Cretan_fans\data'
addpath 'C:\Users\r_ott\Dropbox\Richard\NAGRA\Data\Cosmo'
addpath 'C:\Users\r_ott\Dropbox\Richard\PhD_ETH\matlab\CRONUS cosmo calculation\cronusearth-2.0'
addpath 'C:\Users\r_ott\Dropbox\Richard\PhD_ETH\matlab\InversionBasics\MCMC_book'
addpath '.\subroutines'

% User choice and load data --------------------------------------------- %
scaling_model = 'st';  % scaling model, for nomenclature see CronusCalc
% [num10,txt10,~] = xlsread('10Be_data_CRONUS.xlsx',2);     % load 10Be data
% [num36,txt36,~] = xlsread('36Cl_data_CRONUS.xlsx',2);     % load 36Cl data
[num10,txt10,~] = xlsread('samples.xlsx','10Be Cronus');     % load 10Be data
[num36,txt36,~] = xlsread('samples.xlsx','36Cl Cronus');     % load 36Cl data
[Xdata,~,rawX] = xlsread('samples.xlsx','Samp_comp_for_Matlab_Bedrock');     % load compositional data

soil_mass       = 80;      % average soil mass in g/cm²
DEMdata.method = 'basin';  % Do you want to compute the erosion rate for a specific 'location', or a 'basin'
ind = input('Which of the samples would you like to run (must be the same index in both input tables) ');

% select the desired sample parameters from the input tables
if ind ~= 0
    num10 = num10(ind,:); txt10 = txt10(ind,:);
    num36 = num36(ind,:); txt36 = txt36(ind,:);  % currently use an ind = 2 because this sample is at different positions in my two current input tables
    
    if isnan(rawX{ind+1,2})
        X.fQzS = Xdata(ind,1); 
        X.fCaS = Xdata(ind,2); 
        X.fXS  = Xdata(ind,3);  
        X.mode = 'soil';
    else
        X.fQzB = Xdata(ind,1); 
        X.fCaB = Xdata(ind,2); 
        X.fXB  = Xdata(ind,3);  
        X.mode = 'bedrock';
    end
end
end    

%% assign data and initial basin calculations --------------------------- %

% in case your denudation rate is from an alluvial sample and you desire a
% pixel-by-pixel calculated production rate provide a DEM
if strcmpi('basin',DEMdata.method)
    DEMdata.DEM = GRIDobj();      % interactively choose the DEM that encompasses the basin
    DEMdata.export = 1;           % do you want to save the individual sample scaling factors as .mat file?
    % This can be useful when the computation for scaling schemes like 'sa'
    % and 'sf'  takes a long time for a big basin and you want the scaling
    % factors saved for later
    
    [DEMdata.DB,DEMdata.utmzone] = getBasins(DEMdata.DEM,num10(:,2),num10(:,1),'ll');  % delineate drainage basins and check their geometry
end


%% Calculate production rates ------------------------------------------- %

pars10 = Cronus_prep10(num10,scaling_model,DEMdata);
pars36 = Cronus_prep36(num36,scaling_model,DEMdata);


%% Run MCMC inversion for "real" denudation rate and weathering rate ---- %

% Priors
D = [5,1e3];                           % Denudation min/max in mm/ka
% second prior would be either bedrock or soil quartz content and is
% automatically set to fractions between 0 and 1 in the MCMC function.

% run inversion
[X,MAP,post,W] = pairedCRN_MCMC(pars10,pars36,scaling_model,D,X,soil_mass);


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
disp(['Denudation rate = ' num2str(MAP(2)/pars10.sp10.rb*10) ' mm/ka'])
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
disp(['The calculated weathering rate = ' num2str(W/pars10.sp10.rb*10) ' mm/ka'])

export = input('Do you want to export your results? "y" or "n"? ','s');
if strcmpi(export,'y')
    vars = {'Name','fQzS','fCaS','fXS','fQzB','fCaB','fXB','W','D'};
    out_table = table(txt10,X.fQzS,X.fCaS,X.fXS,X.fQzB,X.fCaB,X.fXB,W/pars10.sp10.rb*10,MAP(2)/pars10.sp10.rb*10 ,'VariableNames',vars);
    writetable(out_table,[ '.\output\JO\pairedN\' txt10{1} '_' X.mode '.xlsx'])
end
