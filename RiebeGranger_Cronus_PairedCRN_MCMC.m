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
[Xdata,~,rawX] = xlsread('samples.xlsx','Sample_composition_for Matlab');     % load compositional data

soil_mass       = 80;      % average soil mass in g/cm²
DEMdata.method = 'basin';  % Do you want to compute the erosion rate for a specific 'location', or a 'basin'
ind = input('Which of the samples would you like to run (must be the same index in both input tables) ');

if ind ~= 0
    num10 = num10(ind,:); txt10 = txt10(ind,:);
    num36 = num36(ind,:); txt36 = txt36(ind,:);  % currently use an ind = 2 because this sample is at different positions in my two current input tables
    
    if isnan(rawX{ind+1,2})
        X.fQzS = Xdata(ind,1); 
        X.fCaS = Xdata(ind,2); 
        X.fXS  = Xdata(ind,3);  
        Xmode = 'soil';
    else
        X.fQzB = Xdata(ind,1); 
        X.fCaB = Xdata(ind,2); 
        X.fXB  = Xdata(ind,3);  
        Xmode = 'bedrock';
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
    
    [DEMdata.DB,DEMdata.utmzone] = getBasins(DEM,num10(:,2),num10(:,1),'ll');  % delineate drainage basins and check their geometry
end

pp=physpars();                               % get physical parameters

%% Calculate production rates ------------------------------------------- %

% First, determine the effective neutron attenuation length following
% Marrero, 2016.
if isnan(num10(13)) || isnan(num36(11))
    Leff = neutron_att_length_DEM(DEM,utmzone);
    num10(13) = Leff;
    num36(11) = Leff;
end

if strcmpi('basin',DEMdata)
    [nominal10,uncerts10,sp10,sf10,cp10,maxage10,maxdepth10,erate_raw10] = Cronus_prep10(num10,...
        scaling_model,pp,DEMdata,DEM,DB,utmzone);
    [nominal36,uncerts36,sp36,sf36,cp36,maxage36,maxdepth36,erate_raw36] = Cronus_prep36(num36,...
        scaling_model,pp,DEMdata,DEM,DB,utmzone);
else
    [nominal10,uncerts10,sp10,sf10,cp10,maxage10,maxdepth10,erate_raw10] = Cronus_prep10(num10,...
        scaling_model,pp,DEMdata);
    [nominal36,uncerts36,sp36,sf36,cp36,maxage36,maxdepth36,erate_raw36] = Cronus_prep36(num36,...
        scaling_model,pp,DEMdata);
end

%% Compute nuclide concentrations for different enrichment/depletions and match measured values ------------- %
tic
% ------------------------------------------------------------------- %
% set some constants
sp10.depthtotop = soil_mass;           % set depth to soil bedrock interface
sp36.depthtotop = soil_mass;           % set depth to soil bedrock interface
soil_depths = 1:0.1:soil_mass; 
% Figure out the maximum possible depth at which we'll ever need a
% production rate.  
maxage=500;                               % 2000ka should be saturated for 36Cl. By trial and error I found that 500 does only ~0.5% diff in the result but doubles the speed
Nobs = [nominal10(9);nominal36(1)];       % measured concentration
dNobs =[uncerts10(9);uncerts36(1)];       %  uncertainty of observation

% set inversion parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Nmax = 1e3;                            % maximum number of models
k = 0.04;                              % universal step size tuned to parameter range
% stds = [0.03;2];                       % step size (crucial parameter)

% PRIORS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
D = [5,1e3];                           % Denudation min/max in mm/ka
D = D/10*sp10.rb;                      % convert to g/cm²/ka for Cronus, I HOPE THIS CONVERSION IS CORRECT
pprior_cur = 0;                        % only flat priors 

switch Xmode
    case 'soil'
        pnames = {'fQzB','D'};                 % names of parameters 
        fQzB = [0,1];                          % fraction of quartz in bedrock   
        range_in = [diff(fQzB);diff(D)];       % ranges of parameters
    case 'bedrock'
        pnames = {'fQzS','D'};                 % names of parameters 
        fQzS = [0,1];                          % fraction of quartz in soil
        range_in = [diff(fQzS);diff(D)];       % ranges of parameters
end


% Resolution, under which parameters resolution does stop the model
err_max = Nobs/100;                    % in at/g for both nuclides, here defined to be 1% of N

nd = length(pnames);                   % number of dimensions

% INIITAL MODEL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
m0 = nan(nd,1);
switch Xmode
    case 'soil'
        m0(1) = fQzB(1)+rand*diff(fQzB);       % create random parameters
    case 'bedrock'
        m0(1) = fQzS(1)+rand*diff(fQzS);       % create random parameters
end
m0(2) = D(1)+rand*diff(D);             % create random parameters

% the covariance matrix should be of the variance (sigma^2) rather than the
% standard deviation (sigma)
covariance_matrix = diag(ones(2,1).*[uncerts10(9);uncerts36(1)].^2);    % set up covariance matrix

% Run initial forward model -----------------------------------------------
Xcur  = X;
switch Xmode
    case 'soil'
        Xcur.fQzB = m0(1); 
        Xcur.fXB  = X.fXS * (m0(1)/X.fQzS);    % the other insoluble mineral should behave like quartz
        Xcur.fCaB = 1 - m0(1) - Xcur.fXB;     % after subtracting the insoluble minerals, the rest should be the soluble mineral
    case 'bedrock'
        Xcur.fQzS = m0(1); 
        Xcur.fXS  = X.fXB * (m0(1)/X.fQzB);    % the other insoluble mineral should behave like quartz
        Xcur.fCaS = 1 - m0(1) - Xcur.fXS;     % after subtracting the insoluble minerals, the rest should be the soluble mineral
end

[N10m,N36m] = paired_N_forward(pp,sp10,sp36,sf10,sf36,cp10,cp36,maxage,scaling_model,soil_mass,m0(2),Xcur);
obs_err     = [N10m,N36m]' - Nobs;                          % observational error
ln_like_cur = (-1/2)*sum((obs_err./dNobs).^2);              % ln likelihood of model
ln_prob_cur = ln_like_cur + log(pprior_cur);                % ln probability model

% MAIN LOOP OF INVERSION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
current = m0;
% post = nan(n,length(current)+1);
nacc = 0;
count = 0;
it = 0; 

con = 1;
while con      % run this while loop until modelled values meet stopping criterion
    it = it+1;        

    % make new random parameters
%     candidate = current + randn(2,1).*stds;
    candidate = current + randn(nd,1).* k .*range_in;
    
    switch Xmode
        case 'soil'
            Xcur.fQzB = candidate(1); 
            Xcur.fXB  = X.fXS * (candidate(1)/X.fQzS);
            Xcur.fCaB = 1 - candidate(1) - Xcur.fXB;        
            % if random parameters are outside of prior range get a new candidate
            while any([candidate < [fQzB(1);D(1)] ; candidate > [fQzB(2);D(2)]; sum(Xcur.fQzB+Xcur.fCaB+Xcur.fXB) > 1.001])
                candidate = current + randn(nd,1).* k .*range_in;
            end
        case 'bedrock'
            Xcur.fQzS = candidate(1); 
            Xcur.fXS  = X.fXB * (candidate(1)/X.fQzB);
            Xcur.fCaS = 1 - candidate(1) - Xcur.fXS;
            while any([candidate < [fQzS(1);D(1)] ; candidate > [fQzS(2);D(2)]; sum(Xcur.fQzS+Xcur.fCaS+Xcur.fXS) > 1.001])
                candidate = current + randn(nd,1).* k .*range_in;
            end
    end

    
    % calculate logpriors for unifrom distribution
    pprior_cand = 0;   % we're using uniform priors

    % new forward model ------------------------------------------------- %

    
    [N10m,N36m,~,~] = paired_N_forward(pp,sp10,sp36,sf10,sf36,cp10,cp36,maxage,scaling_model,soil_mass,candidate(2),Xcur);
    obs_err     = [N10m,N36m]' - Nobs;                          % observational error
    ln_like_cand = (-1/2)*sum((obs_err./dNobs).^2);             % ln likelihood of model    
    lr1 = (-1/2)*sum((candidate-current).^2./(k .*range_in).^2);
    lr2 = (-1/2)*sum((current-candidate).^2./(k .*range_in).^2);
    
    ln_alpha = ln_like_cand + pprior_cand +lr1 - pprior_cur - lr2 - ln_like_cur; % probability candidate, no need to log pprior
    
    if (ln_alpha > 0)
        ln_alpha = 0;
    end

    % Generate a U(0,1) random number and take its logarithm.
    ln_t = log(rand());
    
        % Accept or reject the step.
    if (ln_t < ln_alpha)
        % Accept the step.
        current = candidate;
        ln_like_cur = ln_like_cand;
        pprior_cur = pprior_cand;
        nacc = nacc + 1;
        post(nacc,:) = [candidate; ln_like_cand+pprior_cand];
        disp(it)
        if ln_like_cand+pprior_cand > -1
            k = 0.01;
        elseif nacc > 350
            disp('It seems like the algorithm has problems converging and will terminate now')
            MAP = candidate;
            break
        end
    end
    
    if sum(abs(obs_err) < err_max) == 2   % if both modelled nuclide concentration are close to measured values
        con = 0;                          % stop loop
        MAP = candidate;
    else
        con = 1;                          % if values too far off, keep going
    end
    
end
toc

% Calculate weathering rate --------------------------------------------- %
switch Xmode
    case 'soil'
        X.fQzB = MAP(1);        
        X.fXB  = X.fXS * (MAP(1)/X.fQzS);
        X.fCaB = 1 - MAP(1) - Xcur.fXB;        
    case 'bedrock'
        X.fQzS = MAP(1);        
        X.fXS  = X.fXB * (MAP(1)/X.fQzB);
        X.fCaS = 1 - MAP(1) - Xcur.fXS;
end
% take MAP solution and calculate soil erosion and soil denudation rate
[~,~,Es,EsWs] = paired_N_forward(pp,sp10,sp36,sf10,sf36,cp10,cp36,maxage,scaling_model,soil_mass,MAP(2),Xcur);
W = EsWs - Es;   % weathering rate of carbonate minerals should be the difference between soil erosion and soil denudation rate
W = W * X.fCaS;      % the landscape weathering rate should then be the carbonate weathering rate * the carbonate concentration in the soil       
% NEED TO THINK MORE IF THIS LAST STEP IS CORRECT!

%% OUTPUT RESULTS ------------------------------------------------------- %

% plot the MCMC chains 
figure()
subplot(1,3,1); plot(post(:,1))
switch Xmode
    case 'soil'
        xlabel('iteration'); ylabel(pnames{1}); ylim(fQzB);
    case 'bedrock'
        xlabel('iteration'); ylabel(pnames{1}); ylim(fQzS);
end

subplot(1,3,2); plot(post(:,2))
xlabel('iteration'); ylabel(pnames{2}); ylim(D)
subplot(1,3,3); plot(post(:,3))
xlabel('iteration'); ylabel('log posterior probability')

% Report the values
disp(['Denudation rate = ' num2str(MAP(2)/sp10.rb*10) ' mm/ka'])
switch Xmode
    case 'soil'
        disp(['Fraction of quartz in bedrock fQzB = ' num2str(X.fQzB)])
        disp(['Fraction of X in bedrock fXB = ' num2str(X.fXB)])  
        disp(['Fraction of calcite in bedrock fCaB = ' num2str(X.fCaB)])        
    case 'bedrock'
        disp(['Fraction of quartz in soil fQzS = ' num2str(MAP(1))])
        disp(['Fraction of X in soil fXS = ' num2str(X.fXS)])
        disp(['Fraction of calcite in soil fCaS = ' num2str(X.fCaS)])
end
disp(['The calculated weathering rate = ' num2str(W/sp10.rb*10) ' mm/ka'])

export = input('Do you want to export your results? "y" or "n"? ','s');
if strcmpi(export,'y')
    vars = {'Name','fQzS','fCaS','fXS','fQzB','fCaB','fXB','W','D'};
    out_table = table(txt10,X.fQzS,X.fCaS,X.fXS,X.fQzB,X.fCaB,X.fXB,W/sp10.rb*10,MAP(2)/sp10.rb*10 ,'VariableNames',vars);
    writetable(out_table,[ '.\output\' txt10{1} '.xlsx'])
end
