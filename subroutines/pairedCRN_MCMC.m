function [X,MAP,post,W] = pairedCRN_MCMC(pars10,pars36,scaling_model,D,X)
% Calculates the the "real" landscape denudation rate from a paired nuclide
% measurement and soil or bedrock chemistry data.
% Inversion is perfomred with a MCMC.
% Richard Ott, 2021

v2struct(pars10);
v2struct(pars36);
soil_mass = X.soil_mass;

tic
% ------------------------------------------------------------------- %
% set some constants
pp = physpars();                         % get physical parameters 
sp10.depthtotop = soil_mass;             % set depth to soil bedrock interface
sp36.depthtotop = soil_mass;             % set depth to soil bedrock interface

% Figure out the maximum possible depth at which we'll ever need a
% production rate.  
maxage=500;                               % 2000ka should be saturated for 36Cl. By trial and error I found that 500 does only ~0.5% diff in the result but doubles the speed
Nobs = [nominal10(9);nominal36(1)];       % measured concentration
dNobs =[uncerts10(9);uncerts36(1)];       %  uncertainty of observation

% set inversion parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
k = 0.09;                              % universal step size tuned to parameter range 0,04

% PRIORS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
D = D/10*sp10.rb;                      % convert to g/cm²/ka for Cronus, I HOPE THIS CONVERSION IS CORRECT
pprior_cur = 0;                        % only flat priors 

switch X.mode
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
switch X.mode
    case 'soil'
        m0(1) = fQzB(1)+rand*range_in(1);       % create random parameters
    case 'bedrock'
        m0(1) = fQzS(1)+rand*range_in(1);       % create random parameters
end
spini = sp36; spini.depthtotop = 0;            % set depth to top = 0 for initial erosion rate guess
m0(2) = cl36erateraw(pp,spini,sf36,cp36,scaling_model,0);  % start parameters earch at 36Cl rate (for a carbonate composition this should be closer to the reeal denudation)
% m0(2) = D(1)+rand*diff(D);             % create random parameters

% Run initial forward model -----------------------------------------------
counter =0;
Xcur  = X;
switch X.mode
    case 'soil'
        Xcur.fQzB = m0(1); 
        Xcur.fXB  = X.fXS * (m0(1)/X.fQzS);    % the other insoluble mineral should behave like quartz
        Xcur.fCaB = 1 - m0(1) - Xcur.fXB;     % after subtracting the insoluble minerals, the rest should be the soluble mineral
        while any([m0 < [fQzB(1);D(1)] ; m0 > [fQzB(2);D(2)]; sum(Xcur.fQzB+Xcur.fCaB+Xcur.fXB) > 1.001;any([Xcur.fQzB,Xcur.fXB,Xcur.fCaB]<0) ])
            m0(1) = fQzB(1) + rand(1)*range_in(1);
            Xcur.fQzB = m0(1); 
            Xcur.fXB  = X.fXS * (m0(1)/X.fQzS);
            Xcur.fCaB = 1 - m0(1) - Xcur.fXB; 
        end
    case 'bedrock'
        Xcur.fQzS = m0(1); 
        Xcur.fXS  = X.fXB * (m0(1)/X.fQzB);    % the other insoluble mineral should behave like quartz
        Xcur.fCaS = 1 - m0(1) - Xcur.fXS;     % after subtracting the insoluble minerals, the rest should be the soluble mineral
        while any([m0 < [fQzS(1);D(1)] ; m0 > [fQzS(2);D(2)]; sum(Xcur.fQzS+Xcur.fCaS+Xcur.fXS) > 1.001; any([Xcur.fQzS,Xcur.fXS,Xcur.fCaS]<0)])
            m0(1) = fQzS(1) + rand(1)*range_in(1);
            Xcur.fQzS = m0(1); 
            Xcur.fXS  = X.fXB * (m0(1)/X.fQzB);
            Xcur.fCaS = 1 - m0(1) - Xcur.fXS;
            counter = counter+1;
            if counter > 1e5
                error('cant find a decent sample')
            end
        end
end

[N10m,N36m] = paired_N_forward(pp,sp10,sp36,sf10,sf36,cp10,cp36,maxage,scaling_model,soil_mass,m0(2),Xcur);
obs_err     = [N10m,N36m]' - Nobs;                          % observational error
ln_like_cur = (-1/2)*sum((obs_err./dNobs).^2);              % ln likelihood of model
ln_prob_cur = ln_like_cur + log(pprior_cur);                % ln probability model

% MAIN LOOP OF INVERSION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
current = m0;
nacc = 0;
it = 0; 
counter = 0;
con = 1;
while con      % run this while loop until modelled values meet stopping criterion
    it = it+1;        

    % make new random parameters
    candidate = current + randn(nd,1).* k .*range_in;
    
    switch X.mode
        case 'soil'
            Xcur.fQzB = candidate(1); 
            Xcur.fXB  = X.fXS * (candidate(1)/X.fQzS);
            Xcur.fCaB = 1 - candidate(1) - Xcur.fXB;        
            % if random parameters are outside of prior range get a new candidate
            while any([candidate < [fQzB(1);D(1)] ; candidate > [fQzB(2);D(2)]; sum(Xcur.fQzB+Xcur.fCaB+Xcur.fXB) > 1.001;any([Xcur.fQzB,Xcur.fXB,Xcur.fCaB]<0) ])
                candidate = current + randn(nd,1).* k .*range_in;
                Xcur.fQzB = candidate(1); 
                Xcur.fXB  = X.fXS * (candidate(1)/X.fQzS);
                Xcur.fCaB = 1 - candidate(1) - Xcur.fXB; 
                counter = counter+1;
                if counter > 1e5
                    error('cant find a decent sample')
                end
            end
        case 'bedrock'
            Xcur.fQzS = candidate(1); 
            Xcur.fXS  = X.fXB * (candidate(1)/X.fQzB);
            Xcur.fCaS = 1 - candidate(1) - Xcur.fXS;
            counter = 0;
            while any([candidate < [fQzS(1);D(1)] ; candidate > [fQzS(2);D(2)]; sum(Xcur.fQzS+Xcur.fCaS+Xcur.fXS) > 1.001; any([Xcur.fQzS,Xcur.fXS,Xcur.fCaS]<0)])
                candidate = current + randn(nd,1).* k .*range_in;
                Xcur.fQzS = candidate(1); 
                Xcur.fXS  = X.fXB * (candidate(1)/X.fQzB);
                Xcur.fCaS = 1 - candidate(1) - Xcur.fXS;
                counter = counter+1;
                if counter > 1e5
                    error('cant find a decent sample')
                end
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
    
    counter = 0;
end
toc

% Calculate weathering rate --------------------------------------------- %
switch X.mode
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
[~,~,~,W] = paired_N_forward(pp,sp10,sp36,sf10,sf36,cp10,cp36,maxage,scaling_model,soil_mass,MAP(2),Xcur);

% convert denudation rates back to mm/ka
MAP(2) = MAP(2)/sp10.rb*10;
post(:,1) = post(:,1) ./sp10.rb .*10;
W = W/sp10.rb*10;

%save final composition
X = Xcur;
end

