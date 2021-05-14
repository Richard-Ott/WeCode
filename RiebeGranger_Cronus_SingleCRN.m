% This script calculates the denudation rate from bedrock containing a mix
% of insoluble an soluble minerals.
% Loosely based on Riebe and Granger eqn. 14.
% Conversly to Riebe and Granger, 2014, the nuclide concentrations are not
% calculated with exponentials. Nuclide concentraions are calculated used
% CronusCalc (Marrero, 2016). 
% The script calculates the denudation rate for a single nuclide
% measurement of a soluble or an insoluble target mineral. Despite, the
% nuclide cooncentrations, the bedrock or soil chemistry is needed. This code
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
    load 
else  % load your data files
% addpath 'C:\Users\r_ott\Dropbox\Richard\Crete\Cretan_fans\data'
addpath 'C:\Users\r_ott\Dropbox\Richard\NAGRA\Data\Cosmo'
addpath 'C:\Users\r_ott\Dropbox\Richard\PhD_ETH\matlab\CRONUS cosmo calculation\cronusearth-2.0'
addpath 'C:\Users\r_ott\Dropbox\Richard\NAGRA\Data\Water_CH'
addpath 'C:\Users\r_ott\Dropbox\Richard\PhD_ETH\matlab\InversionBasics\MCMC_book'
addpath '.\subroutines'

% User choice and load data --------------------------------------------- %
nuclide = '10Be';      % '10Be' or '36Cl'
scaling_model = 'st';  % scaling model, for nomenclature see CronusCalc
[num,txt,~]    = xlsread('samples.xlsx','10Be Cronus');                    % load CRN data
[Xdata,~,rawX] = xlsread('samples.xlsx','Sample_composition_for Matlab');  % load compositional data
[Wdata,Wtxt,~] = xlsread('Weathering rates.xlsx');                         % load Weathering data


soil_mass       = 80;         % average soil mass in g/cm²
DEMdata.method  = 'basin';    % Do you want to compute the erosion rate for a specific 'location', or a 'basin'
ind = input('Which of the samples would you like to run ');

if ind ~= 0
    num = num(ind,:); txt = txt(ind,:);
    
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
W = Wdata(1,11)*1e3; Wstd = Wdata(1,12)*1e3;
end
switch nuclide; case '10Be'; n = 1; case '36Cl'; n = 2; end    

%% assign data and initial basin calculations --------------------------- %

% in case your denudation rate is from an alluvial sample and you desire a
% pixel-by-pixel calculated production rate provide a DEM
if strcmpi('basin',DEMdata.method)
    DEMdata.DEM = GRIDobj();        % interactively choose the DEM that encompasses the basin
    DEMdata.export = 1;             % do you want to save the individual sample scaling factors as .mat file?
    % This can be useful when the computation for scaling schemes like 'sa'
    % and 'sf'  takes a long time for a big basin and you want the scaling
    % factors saved for later
    
    [DEMdata.DB,DEMdata.utmzone] = getBasins(DEMdata.DEM,num(:,2),num(:,1),'ll');  % delineate drainage basins and check their geometry
end


%% Calculate production rates ------------------------------------------- %

% Calculate production rates
Cronus_prep = {@Cronus_prep10, @Cronus_prep36};

pars = Cronus_prep{n}(num,scaling_model,DEMdata);

%% Run MCMC inversion for "real" denudation rate ------------------------ %

% PRIORS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
D = [5,1e3];                           % Denudation min/max in mm/ka

% run inversion
[X,MAP,post] = singleCRN_MCMC(pars,scaling_model,D,W,X,soil_mass,n);


%% OUTPUT RESULTS ------------------------------------------------------- %
try rho = pars.sp10.rb;  catch rho = pars.sp36.rb; end

% plot the MCMC chains 
figure()
subplot(1,2,1); plot(post(:,1))
xlabel('iteration'); ylabel('Denudation rate g/cm²/ka'); ylim(D)
subplot(1,2,2); plot(post(:,2))
xlabel('iteration'); ylabel('log posterior probability')

% Report the values
disp(['Denudation rate = ' num2str(MAP/rho*10) ' mm/ka'])
switch X.mode
    case 'soil'
        disp(['Fraction of quartz in bedrock fQzB = ' num2str(X.fQzB)])
        disp(['Fraction of X in bedrock fXB = ' num2str(X.fXB)])  
        disp(['Fraction of calcite in bedrock fCaB = ' num2str(X.fCaB)])        
    case 'bedrock'
        disp(['Fraction of quartz in soil fQzS = ' num2str(X.fQzS)])
        disp(['Fraction of X in soil fXS = ' num2str(X.fXS)])
        disp(['Fraction of calcite in soil fCaS = ' num2str(X.fCaS)])
end

export = input('Do you want to export your results? "y" or "n"? ','s');
if strcmpi(export,'y')
    vars = {'Name','fQzS','fCaS','fXS','fQzB','fCaB','fXB','W','D'};
    out_table = table(txt,X.fQzS,X.fCaS,X.fXS,X.fQzB,X.fCaB,X.fXB,W/rho*10,MAP/rho*10 ,'VariableNames',vars);
    writetable(out_table,[ '.\output\JO\10Be\' txt{1} '.xlsx'])
end
