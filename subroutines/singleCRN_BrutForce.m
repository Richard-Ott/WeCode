function [X,MAP,post] = singleCRN_BrutForce(pars,D,X,err)
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
        

% PRIORS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
D = D/10*sp.rb;                        % convert to g/cm²/ka for Cronus
W = W/10*sp.rb;                        % convert to g/cm²/ka 

pnames = {'D'};                        % prior names
range_in = diff(D);                    % ranges of parameters

% Resolution under which parameters resolution does stop the model
err_max = Nobs*err/100;                % in at/g for nuclides
nd = length(pnames);                   % number of dimensions

Drates = linspace(D(1),D(2));
Nm = nan(100,1);
forward_model = {@N10_forward,@N36_forward};        % to avoid opening more switch statements, I compile the functions for each case here
for i = 1:100   

    m0 = Drates(i);
    Xcur = X;
    fE = 1 - W/m0;                          % fraciotn of erosion (to total denudation)
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

    % new forward model ------------------------------------------------- %
    [Nm(i),~] = forward_model{n}(pp,sp,sf,cp,maxage,scaling_model,soil_mass,m0,Xcur);
    obs_err     = Nm - Nobs;                                    % observational error


end

plot(Drates./sp.rb*10,Nm);
hold on
yline(Nobs);


%save final composition
X = Xcur;

% convert denudation rates from g/cm2/ka to mm/ka
MAP = MAP/sp.rb*10;   
post(:,1) = post(:,1) ./sp.rb .*10;

end

