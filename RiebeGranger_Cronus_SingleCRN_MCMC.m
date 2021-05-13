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


soil_mass       = 80;      % average soil mass in g/cm²
DEMdata = 'basin';         % Do you want to compute the erosion rate for a specific 'location', or a 'basin'
ind = input('Which of the samples would you like to run ');

if ind ~= 0
    num = num(ind,:); txt = txt(ind,:);
    
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
W = Wdata(1,11)*1e3; Wstd = Wdata(1,12)*1e3;
switch nuclide; case '10Be'; n = 1; case '36Cl'; n = 2; end    

%% assign data and initial basin calculations --------------------------- %

% in case your denudation rate is from an alluvial sample and you desire a
% pixel-by-pixel calculated production rate provide a DEM
if strcmpi('basin',DEMdata)
    DEM = GRIDobj();        % interactively choose the DEM that encompasses the basin
    export = 1;             % do you want to save the individual sample scaling factors as .mat file?
    % This can be useful when the computation for scaling schemes like 'sa'
    % and 'sf'  takes a long time for a big basin and you want the scaling
    % factors saved for later
    
    [DB,utmzone] = getBasins(DEM,num(:,2),num(:,1),'ll');  % delineate drainage basins and check their geometry
end

pp=physpars();                               % get physical parameters
nsamples=size(num,1);                        % number of samples

%% Calculate production rates ------------------------------------------- %

% First, determine the effective neutron attenuation length following
% Marrero, 2016.
if isnan(num(13)) && strcmpi(nuclide,'10Be')
    Leff = neutron_att_length_DEM(DEM,utmzone);
    num(13) = Leff;
elseif isnan(num(11)) && strcmpi(nuclide,'36Cl')
    Leff = neutron_att_length_DEM(DEM,utmzone);
    num(11) = Leff;
end

% Calculate production rates
Cronus_prep = {@Cronus_prep10, @Cronus_prep36};
if strcmpi('basin',DEMdata)
    [nominal,uncerts,sp,sf,cp,maxage,maxdepth,erate_raw] = Cronus_prep{n}(num,...
        scaling_model,pp,DEMdata,DEM,DB,utmzone);
else
    [nominal,uncerts,sp,sf,cp,maxage,maxdepth,erate_raw] = Cronus_prep{n}(num,...
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
maxage=500;                            % 2000ka should be saturated for 36Cl. By trial and error I found that 500 does only ~0.5% diff in the result but doubles the speed
data_ind = n+8-9*(n-1);                % index I use for referencing into the nominal and uncerts arrays, because unfortunately Cronus uses different order for the data depending on the nuclide
Nobs = nominal(data_ind);              % measured concentration
dNobs =uncerts(data_ind);              % uncertainty of observation
        
% set inversion parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Nmax = 1e3;                            % maximum number of models
k = 0.04;                              % universal step size tuned to parameter range

% PRIORS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
D = [5,1e3];                           % Denudation min/max in mm/ka
D = D/10*sp.rb;                        % convert to g/cm²/ka for Cronus, I HOPE THIS CONVERSION IS CORRECT
W = W/10*sp.rb;                        % convert to g/cm²/ka 
pprior_cur = 0;                        % only flat priors 

pnames = {'D'};
range_in = diff(D);                    % ranges of parameters

% Resolution, under which parameters resolution does stop the model
err_max = Nobs/100;                    % in at/g for nuclides, here defined to be 1% of N
nd = length(pnames);                   % number of dimensions

% INIITAL MODEL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
m0 = D(1)+rand*diff(D);                % create random parameters

% the covariance matrix should be of the variance (sigma^2) rather than the
% standard deviation (sigma)
covariance_matrix = diag(ones(1,1).*uncerts(data_ind).^2);    % set up covariance matrix

% Run initial forward model -----------------------------------------------
Xcur  = X;
fE = 1 - W/m0;                          % fraciotn of erosion (to total denudation)
switch Xmode
    case 'soil'
        R  = (Xcur.fQzS + Xcur.fXS)/Xcur.fCaS;     % ratio of insoluble to soluble material
        Xcur.fQzB = (R*fE) / (1+ R*fE + R*fE*(X.fXS/X.fQzS) + X.fXS/X.fQzS); 
        Xcur.fXB  = Xcur.fQzB* (X.fXS/X.fQzS);  
        Xcur.fCaB = 1 - Xcur.fQzB - Xcur.fXB;
    case 'bedrock'
        R  = (Xcur.fQzB + Xcur.fXB)/Xcur.fCaB;     % ratio of insoluble to soluble material
        Xcur.fQzS = (R/fE) / (1+ R/fE + R/fE*(X.fXB/X.fQzB) + X.fXB/X.fQzB); 
        Xcur.fXS  = Xcur.fQzS* (X.fXB/X.fQzB);  
        Xcur.fCaS = 1 - Xcur.fQzS - Xcur.fXS;
end

forward_model = {@N10_forward,@N36forward};        % to avoid opening more switch statements, I compile the functions for each case here

[Nm,~] = forward_model{n}(pp,sp,sf,cp,maxage,scaling_model,soil_mass,m0,Xcur);
obs_err     = Nm - Nobs;                           % observational error
ln_like_cur = (-1/2)*sum((obs_err./dNobs).^2);     % ln likelihood of model
ln_prob_cur = ln_like_cur + log(pprior_cur);       % ln probability model

% MAIN LOOP OF INVERSION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
current = m0;
nacc = 0;
count = 0;
it = 0; 

con = 1;
while con      % run this while loop until modelled values meet stopping criterion
    it = it+1;        

    % make new random parameters
    candidate = current + randn(nd,1).* k .*range_in;
    
    fE = 1 - W/candidate;                          % fraciotn of erosion (to total denudation)
    switch Xmode
        case 'soil'
            R  = (Xcur.fQzS + Xcur.fXS)/Xcur.fCaS;     % ratio of insoluble to soluble material
            Xcur.fQzB = (R*fE) / (1+ R*fE + R*fE*(X.fXS/X.fQzS) + X.fXS/X.fQzS); 
            Xcur.fXB  = Xcur.fQzB* (X.fXS/X.fQzS);  
            Xcur.fCaB = 1 - Xcur.fQzB - Xcur.fXB;  
            % if random parameters are outside of prior range get a new candidate
            while any([candidate < D(1) ; candidate > D(2)])
                candidate = current + randn(nd,1).* k .*range_in;
            end
        case 'bedrock'
            R  = (Xcur.fQzB + Xcur.fXB)/Xcur.fCaB;     % ratio of insoluble to soluble material
            Xcur.fQzS = (R/fE) / (1+ R/fE + R/fE*(X.fXB/X.fQzB) + X.fXB/X.fQzB); 
            Xcur.fXS  = Xcur.fQzS* (X.fXB/X.fQzB);  
            Xcur.fCaS = 1 - Xcur.fQzS - Xcur.fXS;
            while any([candidate < D(1); candidate > D(2)])
                candidate = current + randn(nd,1).* k .*range_in;
            end
    end
    

    % calculate logpriors for unifrom distribution
    pprior_cand = 0;   % we're using uniform priors

    % new forward model ------------------------------------------------- %
    [Nm,~] = forward_model{n}(pp,sp,sf,cp,maxage,scaling_model,soil_mass,candidate,Xcur);
    obs_err     = Nm - Nobs;                                    % observational error
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
    
    if abs(obs_err) < err_max   % if the modelled nuclide concentration is close to measured values
        con = 0;                % stop loop
        MAP = candidate;
    else
        con = 1;                % if values too far off, keep going
    end
    
end
toc

X = Xcur;

%% OUTPUT RESULTS ------------------------------------------------------- %

% plot the MCMC chains 
figure()
subplot(1,2,1); plot(post(:,1))
xlabel('iteration'); ylabel(pnames); ylim(D)
subplot(1,2,2); plot(post(:,2))
xlabel('iteration'); ylabel('log posterior probability')

% Report the values
disp(['Denudation rate = ' num2str(MAP/sp.rb*10) ' mm/ka'])
switch Xmode
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
    out_table = table(txt10,X.fQzS,X.fCaS,X.fXS,X.fQzB,X.fCaB,X.fXB,W/sp.rb*10,MAP/sp.rb*10 ,'VariableNames',vars);
    writetable(out_table,[ '.\output\' txt10{1} '.xlsx'])
end
