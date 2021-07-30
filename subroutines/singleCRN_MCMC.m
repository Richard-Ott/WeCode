function [X,MAP,post] = singleCRN_MCMC(pars,scaling_model,D,X,err)
% Calculates the "real" denudation rate from a nuclide measurement, the
% bedrock or soil chemistry and a weathering rate.
% The solution is found through a MCMC algorithm.
% Richard Ott, 2021

v2struct(pars)   % unpack variables in parameters structure
W = X.W;
Wstd = X.Wstd;
soil_mass = X.soil_mass;
n = X.n;

% I have to put this ugly if-else-statement here to rename the variables
% otherwise I'd need to define new parameter calculation funtions...
if exist('sp10')
    nominal = nominal10; uncerts = uncerts10; sp = sp10; sf = sf10; cp = cp10;
    clear nominal10 uncerts10 sp10 sf10 cp10
elseif exist('sp36')
    nominal = nominal36; uncerts = uncerts36; sp = sp36; sf = sf36; cp = cp36;
    clear nominal36 uncerts36 sp36 sf36 cp36
end
    
tic
% ------------------------------------------------------------------- %
% set some constants
pp = physpars();
sp.depthtotop = soil_mass;           % set depth to soil bedrock interface
% Figure out the maximum possible depth at which we'll ever need a
% production rate.  
maxage=500;                            % 2000ka should be saturated for 36Cl. By trial and error I found that 500 does only ~0.5% diff in the result but doubles the speed
data_ind = n+8-9*(n-1);                % index I use for referencing into the nominal and uncerts arrays, because unfortunately Cronus uses different order for the data depending on the nuclide
Nobs = nominal(data_ind);              % measured concentration
dNobs =uncerts(data_ind);              % uncertainty of observation
        
% set inversion parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
k = 0.04;                              % universal step size tuned to parameter range

% PRIORS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
D = D/10*sp.rb;                        % convert to g/cm²/ka for Cronus, I HOPE THIS CONVERSION IS CORRECT
W = W/10*sp.rb;                        % convert to g/cm²/ka 
Wstd = Wstd/10*sp.rb;                  % convert to g/cm²/ka 
pprior_cur = 0;                        % only flat priors 

pnames = {'D'};
range_in = diff(D);                    % ranges of parameters

% Resolution, under which parameters resolution does stop the model
err_max = Nobs*err/100;                    % in at/g for nuclides, here defined to be 1% of N
nd = length(pnames);                   % number of dimensions

% INIITAL MODEL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
m0 = D(1)+rand*diff(D);                % create random parameters

% the covariance matrix should be of the variance (sigma^2) rather than the
% standard deviation (sigma)
covariance_matrix = diag(ones(1,1).*uncerts(data_ind).^2);    % set up covariance matrix

% Run initial forward model -----------------------------------------------
Xcur  = X;
fE = 1 - W/m0;                          % fraciotn of erosion (to total denudation)
switch X.mode
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

forward_model = {@N10_forward,@N36_forward};        % to avoid opening more switch statements, I compile the functions for each case here

[Nm,~] = forward_model{n}(pp,sp,sf,cp,maxage,scaling_model,soil_mass,m0,Xcur);
obs_err     = Nm - Nobs;                           % observational error
ln_like_cur = (-1/2)*sum((obs_err./dNobs).^2);     % ln likelihood of model
ln_prob_cur = ln_like_cur + log(pprior_cur);       % ln probability model

% MAIN LOOP OF INVERSION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
current = m0;
nacc = 0;
it = 0; 

con = 1;
up = 0;
while con      % run this while loop until modelled values meet stopping criterion
    it = it+1;        

    % make new random parameters ------------------------------------------
    candidate = current + randn(nd,1).* k .*range_in;
    while up == 1 && (candidate < current)   % if the last accepted denudation rate resulted in too high concentration we have to further increase denudation
        candidate = current + randn(nd,1).* k .*range_in;   % this loops saves us from calculting some useless forward models
    end
    
    fE = 1 - W/candidate;                          % fraciotn of erosion (to total denudation)
    switch X.mode
        case 'soil'
            R  = (Xcur.fQzS + Xcur.fXS)/Xcur.fCaS;             % ratio of insoluble to soluble material
            Xcur.fQzB = (R*fE) / (1+ R*fE + R*fE*(X.fXS/X.fQzS) + X.fXS/X.fQzS); 
            Xcur.fXB  = Xcur.fQzB* (X.fXS/X.fQzS);  
            Xcur.fCaB = 1 - Xcur.fQzB - Xcur.fXB;  
            % if random parameters are outside of prior range get a new candidate
            while any([candidate < D(1) ; candidate > D(2); sum(Xcur.fQzB+Xcur.fCaB+Xcur.fXB) > 1.001;any([Xcur.fQzB,Xcur.fXB,Xcur.fCaB]<0) ])
                candidate = current + randn(nd,1).* k .*range_in;
                fE = 1 - W/candidate;                          % fraciotn of erosion (to total denudation)
                R  = (Xcur.fQzS + Xcur.fXS)/Xcur.fCaS;         % ratio of insoluble to soluble material
                Xcur.fQzB = (R*fE) / (1+ R*fE + R*fE*(X.fXS/X.fQzS) + X.fXS/X.fQzS); 
                Xcur.fXB  = Xcur.fQzB* (X.fXS/X.fQzS);  
                Xcur.fCaB = 1 - Xcur.fQzB - Xcur.fXB;  
            end
        case 'bedrock'
            R  = (Xcur.fQzB + Xcur.fXB)/Xcur.fCaB;             % ratio of insoluble to soluble material
            Xcur.fQzS = (R/fE) / (1+ R/fE + R/fE*(X.fXB/X.fQzB) + X.fXB/X.fQzB); 
            Xcur.fXS  = Xcur.fQzS* (X.fXB/X.fQzB);  
            Xcur.fCaS = 1 - Xcur.fQzS - Xcur.fXS;
            while any([candidate < D(1); candidate > D(2); sum(Xcur.fQzS+Xcur.fCaS+Xcur.fXS) > 1.001;any([Xcur.fQzS,Xcur.fXS,Xcur.fCaS]<0)])
                candidate = current + randn(nd,1).* k .*range_in;
                fE = 1 - W/candidate;                          % fraciotn of erosion (to total denudation)
                R  = (Xcur.fQzB + Xcur.fXB)/Xcur.fCaB;         % ratio of insoluble to soluble material
                Xcur.fQzS = (R/fE) / (1+ R/fE + R/fE*(X.fXB/X.fQzB) + X.fXB/X.fQzB); 
                Xcur.fXS  = Xcur.fQzS* (X.fXB/X.fQzB);  
                Xcur.fCaS = 1 - Xcur.fQzS - Xcur.fXS;
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

    % Generate a U(0,1) random number and take its logarithm 
    ln_t = log(rand());
    
        % Accept or reject the step ---------------------------------------
    if (ln_t < ln_alpha)
        % Accept the step.
        current = candidate;
        ln_like_cur = ln_like_cand;
        pprior_cur = pprior_cand;
        nacc = nacc + 1;
        post(nacc,:) = [candidate; ln_like_cand+pprior_cand];
        if obs_err > 0; up = 1; else up = 0; end   % we know if the modelled N was to high we have to increase denudation rate
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

%save final composition
X = Xcur;

% convert denudation rates back to mm/ka
MAP = MAP/sp.rb*10;
post(:,1) = post(:,1) ./sp.rb .*10;

end

