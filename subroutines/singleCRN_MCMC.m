function [X,MAP,post] = singleCRN_MCMC(pars,D,X,err)
% (not recommended, use singleCRN_Optim instead, it is much faster)
% Calculates the "real" denudation rate from a nuclide measurement, the
% bedrock or soil chemistry and a weathering rate.
% The solution is found through a MCMC algorithm.
% Richard Ott, 2021

global scaling_model
v2struct(pars)   % unpack variables in parameters structure
W = X.W;
soil_mass = X.soil_mass;
n = X.n;

% I have to put this ugly if-else-statement here to rename the variables
% otherwise I'd need to define new parameter calculation funtions...
if exist('sp10')
    nominal = nominal10; uncerts = uncerts10; sp = sp10; sf = sf10; cp = cp10; maxage = maxage10; 
    clear nominal10 uncerts10 sp10 sf10 cp10 maxage10
elseif exist('sp36')
    nominal = nominal36; uncerts = uncerts36; sp = sp36; sf = sf36; cp = cp36; maxage = maxage36; 
    clear nominal36 uncerts36 sp36 sf36 cp36 maxage36
end
    
tic
% ------------------------------------------------------------------- %
% set some constants
pp = physpars();
sp.depthtotop = soil_mass;           % set depth to soil bedrock interface
% Figure out the maximum possible depth at which we'll ever need a
% production rate.  
% maxage=500;                          % 2000ka should be saturated for 36Cl. By trial and error I found that 500 does only ~0.5% diff in the result but doubles the speed
data_ind = n+8-9*(n-1);                % index I use for referencing into the nominal and uncerts arrays, because unfortunately Cronus uses different order for the data depending on the nuclide
Nobs = nominal(data_ind);              % measured concentration
dNobs =uncerts(data_ind);              % uncertainty of observation

% Minimum denudation rate %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The enrichment cannot produce denudation rates that would lead to
% enrichment of insoluble minerals that go beyonf 100% of the soil fraction
switch X.mode; case 'soil'; insol = X.fQzS + X.fXS; case 'bedrock'; insol = X.fQzB + X.fXB; end
if insol > (D(1)-W)/D(1); D(1) = W/(1-insol); end
        
% set inversion parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
k = 0.08;                              % universal step size tuned to parameter range, larger number -> larger steps
% This step size determines how "jumpy" the random walk is and might need
% to be tuned if the alogirthm is not working well for your samples

% PRIORS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
D = D/10*sp.rb;                        % convert to g/cm²/ka for Cronus
W = W/10*sp.rb;                        % convert to g/cm²/ka 
ln_pprior_cur = 0;                        % only flat priors 

pnames = {'D'};                        % prior names
range_in = diff(D);                    % ranges of parameters

% Resolution under which parameters resolution does stop the model
err_max = Nobs*err/100;                % in at/g for nuclides
nd = length(pnames);                   % number of dimensions

% INIITAL MODEL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% You can either start the parameter search at a random location (first
% line), or at the conventional erosion rate (2nd-4th line). To
% change the method just uncomment and comment the appropriate lines.
% m0 = D(1)+rand*diff(D);            
spini = sp; spini.depthtotop = 0;            % set depth to top = 0 for initial erosion rate guess
Dini = {@be10erateraw, @cl36erateraw};
m0 = Dini{n}(pp,spini,sf,cp,scaling_model,0);


% Run initial forward model -----------------------------------------------
Xcur  = X;
fE = 1 - W/m0;                          % fraction of erosion (to total denudation)
switch X.mode
    case 'soil'
        Xcur.fQzB = Xcur.fQzS * fE;           % quartz fraction 
        Xcur.fXB  = Xcur.fXS  * fE;           % clay fraction
        Xcur.fCaB = 1 - Xcur.fQzB - Xcur.fXB;
    case 'bedrock'
        Xcur.fQzS = Xcur.fQzB * (1/fE);
        Xcur.fXS  = Xcur.fXB  * (1/fE); 
        Xcur.fCaS = 1 - Xcur.fQzS - Xcur.fXS;
end

forward_model = {@N10_forward,@N36_forward};        % to avoid opening more switch statements, I compile the functions for each case here

[Nm,~] = forward_model{n}(pp,sp,sf,cp,maxage,scaling_model,soil_mass,m0,Xcur);
obs_err     = Nm - Nobs;                           % observational error
% if obs_err > 0; up = true; else up = false;end   % OPTIONAL (if uncomment also uncomment the other OPTIONAL lines: check if initial denudation rate is too high or too low
ln_like_cur = (-1/2)*sum((obs_err./dNobs).^2);     % ln likelihood of model
ln_prob_cur = ln_like_cur + ln_pprior_cur;       % ln probability model

% MAIN LOOP OF INVERSION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
current = m0;
nacc = 0;
it = 0; 

con = 1;
counter = 0;
while con      % run this while loop until modelled values meet stopping criterion
    it = it+1;        

    % sample new random parameter ------------------------------------------
    candidate = truncnormrnd(1,current,k .*range_in, D(1),D(2));
    % the if-selse-statement makes sure that the new sample is either
    % larger or smaller than the previous one depending if the last
    % predicted nuclide concentration was too low or too high The lines below 
    % make the code run faster but can lead to problems if you're close to 
    % pure weathering.
%     if up
%         candidate = truncnormrnd(1,current,k .*range_in, current,D(2));
%     else
%         candidate = truncnormrnd(1,current,k .*range_in, D(1),current);
%     end
    
    fE = 1 - W/candidate;                          % fraciotn of erosion (to total denudation)
    switch X.mode
        case 'soil'
            Xcur.fQzB = Xcur.fQzS * fE;           % quartz fraction 
            Xcur.fXB  = Xcur.fXS  * fE;           % clay fraction
            Xcur.fCaB = 1 - Xcur.fQzB - Xcur.fXB;
        case 'bedrock'
            Xcur.fQzS = Xcur.fQzB * (1/fE);
            Xcur.fXS  = Xcur.fXB  * (1/fE); 
            Xcur.fCaS = 1 - Xcur.fQzS - Xcur.fXS;
    end
    
    % calculate logpriors for unifrom distribution
    ln_pprior_cand = 0;   % we're using uniform priors

    % new forward model ------------------------------------------------- %
    [Nm,~] = forward_model{n}(pp,sp,sf,cp,maxage,scaling_model,soil_mass,candidate,Xcur);
    obs_err     = Nm - Nobs;                                    % observational error
    ln_like_cand = (-1/2)*sum((obs_err./dNobs).^2);             % ln likelihood of model    
    lr1 = (-1/2)*sum((candidate-current).^2./(k .*range_in).^2);
    lr2 = (-1/2)*sum((current-candidate).^2./(k .*range_in).^2);
    
    ln_alpha = ln_like_cand + ln_pprior_cand +lr1 - ln_pprior_cur - lr2 - ln_like_cur; % probability candidate, no need to log pprior
    
    if (ln_alpha > 0)
        ln_alpha = 0;
    end

    % Generate a U(0,1) random number and take its logarithm 
    ln_t = log(rand());
    
        % Accept or reject the step ---------------------------------------
    if (ln_t < ln_alpha)
        % Accept the step.
        current = candidate;
        ln_like_cur = ln_like_cand;
        pprior_cur = ln_pprior_cand;
        nacc = nacc + 1;
        post(nacc,:) = [candidate; ln_like_cand+ln_pprior_cand];
%         if obs_err > 0; up = true; else up = false; end   % OPTIONAL: we know if the modelled N was to high we have to increase denudation rate
        disp(it)
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % The k-parameters below might need to be tuned to your samples for
        % efficiency
        if ln_like_cand+ln_pprior_cand > -0.1        % decrease step size once we're close to a solution
            k = 0.003;
        elseif ln_like_cand+ln_pprior_cand > -1
            k = 0.01;
        elseif ln_like_cand+ln_pprior_cand > -2
            k = 0.025;
        elseif nacc > 350
            MAP = candidate;
            error('It seems like the algorithm has problems converging and will terminate now. Rerun the algorithm, check the input, or play with the inversion parameters (e.g., k)')
        end
    end
    
    if abs(obs_err) < err_max   % if the modelled nuclide concentration is close to measured values
        con = 0;                % stop loop
        MAP = candidate;
    else
        con = 1;                % if values too far off, keep going
    end
    counter = 0;
    
end
toc

%save final composition
X = Xcur;

% convert denudation rates from g/cm2/ka to mm/ka
MAP = MAP/sp.rb*10;   
post(:,1) = post(:,1) ./sp.rb .*10;

end

