% This script calculates the denudation rate from bedrock containing a mix
% of insoluble an soluble minerals.
% Loosely based on Riebe and Granger eqn. 14.
% Conversly to Riebe and Granger, 2014, the nuclide concentrations are not
% calculated with exponentials. Nuclide concentraions are calculated used
% CronusCalc (Marrero, 2016). 
% The script calculates the denudation rate for a paired nuclide
% measurement of a soluble and an insoluble target mineral. Despite, the
% nuclide cooncentrations, the bedrock chemistry is needed. This code
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
addpath 'C:\Users\r_ott\Dropbox\Richard\Crete\Cretan_fans\data'
addpath 'C:\Users\r_ott\Dropbox\Richard\PhD_ETH\matlab\CRONUS cosmo calculation\cronusearth-2.0'
addpath 'C:\Users\r_ott\Dropbox\Richard\PhD_ETH\matlab\InversionBasics\MCMC_book'
addpath '.\subroutines'

% User choice ----------------------------------------------------------- %
scaling_model = 'st';  % scaling model, for nomenclature see CronusCalc
[num10,txt10,~] = xlsread('10Be_data_CRONUS.xlsx',2);     % load 10Be data
[num36,txt36,~] = xlsread('36Cl_data_CRONUS.xlsx',2);     % load 36Cl data
soil_mass       = 50;                                     % average soil mass in g/cm²

DEMdata = 'location';     % Do you want to compute the erosion rate for a specific 'location', or a 'basin'
ind = input('Which of the samples would you like to run (must be the same index in both input tables?) ');
if ind ~= 0
    num10 = num10(ind,:); txt10 = txt10(ind,:);
    num36 = num36(ind+2,:); txt36 = txt36(ind+2,:);  % currently use an ind = 2
end

X.fQzB = num10(:,35); 
X.fCaB = num36(:,82); 
X.fXB  = num10(:,37);  


%% assign data and initial basin calculations --------------------------- %

% in case your denudation rate is from an alluvial sample and you desire a
% pixel-by-pixel calculated production rate provide a DEM
if strcmpi('basin',DEMdata)
    DEM = GRIDobj();        % interactively choose the DEM that encompasses the basin
    export = 1;             % do you want to save the individual sample scaling factors as .mat file?
    % This can be useful when the computation for scaling schemes like 'sa'
    % and 'sf'  takes a long time for a big basin and you want the scaling
    % factors saved for later
    
    [DB,utmzone] = getBasins(DEM,num10(:,2),num10(:,1),'ll');  % delineate drainage basins and check their geometry
end

pp=physpars();                               % get physical parameters
nsamples=size(num10,1);                      % number of samples

%% Calculate production rates ------------------------------------------- %

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
Nobs = [nominal10(9);nominal36(1)];                       % measured concentration
dNobs =[uncerts10(9);uncerts36(1)];      %  uncertainty of observation

% set inversion parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Nmax = 1e3;                            % maximum number of models
k = 0.04;                              % universal step size tuned to parameter range
% stds = [0.03;2];                       % step size (crucial parameter)

% PRIORS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pnames = {'fQzS','D'};
fQzS = [0,1];                          % fraction of quartz in soil
D = [5,1e3];                           % Denudation min/max in mm/ka
D = D/10*sp10.rb;                      % convert to g/cm²/ka for Cronus, I HOPE THIS CONVERSION IS CORRECT
pprior_cur = 0;                        % only flat priors 
range_in = [diff(fQzS);diff(D)];       % ranges of parameters

% Resolution, under which parameters resolution does stop the model
err_max = Nobs/100;                    % in at/g for both nuclides, here defined to be 1% of N

nd = length(pnames);                   % number of dimensions

% INIITAL MODEL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
m0 = nan(nd,1);
m0(1) = fQzS(1)+rand*diff(fQzS);       % create random parameters
m0(2) = D(1)+rand*diff(D);             % create random parameters

% the covariance matrix should be of the variance (sigma^2) rather than the
% standard deviation (sigma)
covariance_matrix = diag(ones(2,1).*[uncerts10(9);uncerts36(1)].^2);    % set up covariance matrix

% Run initial forward model -----------------------------------------------
Xcur  = X;
Xcur.fQzS = m0(1); 
Xcur.fXS  = X.fXB * (m0(1)/X.fQzB);    % the other insoluble mineral should behave like quartz
Xcur.fCaS = 1 - m0(1) - Xcur.fXS;     % after subtracting the insoluble minerals, the rest should be the soluble mineral
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
    
    % if random parameters are outside of prior range get a new candidate
    while any(candidate > [fQzS(1);D(1)] & candidate < [fQzS(1);D(1)])
        candidate = current + randn(nd,1).* k .*range_in;
    end
        
    
    % calculate logpriors for unifrom distribution
    pprior_cand = 0;   % we're using uniform priors

    % new forward model ------------------------------------------------- %
    Xcur.fQzS = candidate(1); 
    Xcur.fXS  = X.fXB * (candidate(1)/X.fQzB);
    Xcur.fCaS = 1 - candidate(1) - Xcur.fXS;
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
Xcur.fQzS = MAP(1);        
Xcur.fXS  = X.fXB * (MAP(1)/X.fQzB);
Xcur.fCaS = 1 - MAP(1) - Xcur.fXS;
% take MAP solution and calculate soil erosion and soil denudation rate
[~,~,Es,EsWs] = paired_N_forward(pp,sp10,sp36,sf10,sf36,cp10,cp36,maxage,scaling_model,soil_mass,MAP(2),Xcur);
W = EsWs - Es;   % weathering should be the difference between soil erosion and soil denudation rate

%% OUTPUT RESULTS ------------------------------------------------------- %

% plot the chains
figure()
subplot(1,3,1); plot(post(:,1))
xlabel('iteration'); ylabel(pnames{1}); ylim(fQzS)
subplot(1,3,2); plot(post(:,2))
xlabel('iteration'); ylabel(pnames{2}); ylim(D)
subplot(1,3,3); plot(post(:,3))
xlabel('iteration'); ylabel('log posterior probability')

% Report the values
disp(['Denudation rate = ' num2str(MAP(2)/sp10.rb*10) ' mm/ka'])
disp(['Fraction of quartz in soil fQzS = ' num2str(MAP(1))])
disp(['Fraction of calcite in soil fCaS = ' num2str(1-MAP(1))])
disp(['The calculated weathering rate = ' num2str(W/sp10.rb*10) ' mm/ka'])


