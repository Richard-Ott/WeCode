function [D_Nmax,Nmax] = solCRN_D_Nmax(pars,D,X,err)
% Finds the maximum nuclide concentration and its denudation rate for a
% soluble target mineral.
% Input: 
%       - pars: Cronus parameters for CRN calculation see (Cronus_prepXX.m)
%       - D: [Dmin, Dmax] defines search boundaries for denudation rate
%           in mm/ka
%       - X: compositional structure (see CosmoDataRead.m)
%       - err: threshold of nuclide concentration error at which 
%           optimization has converged in % of N
% Richard Ott, 2021

global scaling_model
v2struct(pars)   % unpack variables in parameters structure
W = X.W;
soil_mass = X.soil_mass;
n = X.n;

% I have to put this ugly if-else-statement here to rename the variables
% otherwise I'd need to define new parameter calculation funtions...
if exist('sp10')
    nominal = nominal10; sp = sp10; sf = sf10; cp = cp10; maxage = maxage10; 
    clear nominal10 uncerts10 sp10 sf10 cp10 maxage10
elseif exist('sp36')
    nominal = nominal36; sp = sp36; sf = sf36; cp = cp36; maxage = maxage36; 
    clear nominal36 uncerts36 sp36 sf36 cp36 maxage36
end
    
tic
% ------------------------------------------------------------------- %
% set some constants
pp = physpars();
sp.depthtotop = soil_mass;             % set depth to soil bedrock interface
data_ind = n+8-9*(n-1);                % index I use for referencing into the nominal and uncerts arrays, because unfortunately Cronus uses different order for the data depending on the nuclide
Nobs = nominal(data_ind);              % measured concentration

D = D/10*sp.rb;                        % convert to g/cm²/ka for Cronus
W = W/10*sp.rb;                        % convert to g/cm²/ka 
forward_model = {@N10_forward,@N36_forward};% to avoid opening more switch statements, I compile the functions for each case here


%% Minimum denudation rate %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The enrichment cannot produce denudation rates that would lead to
% enrichment of insoluble minerals that go beyond 100% of the soil fraction
switch X.mode; case 'soil'; insol = X.fQzS + X.fXS; case 'bedrock'; insol = X.fQzB + X.fXB; end
if insol > (D(1)-W)/D(1)
    D(1) = W/(1-insol);    % mimimum denudation rate Dmin
end

%% Identify maximum possible nuclide concentration for this sample and model
tolerance = Nobs*err/100;            % in at/g for nuclides
% Define options for optimization  
options = optimset('MaxIter',5e4,'TolFun',tolerance,'TolX',0.05); 

x0 = D(1)+3;       % start at denudation rate slightly higher than minimum
fun = @(x) (-1)*Comp_and_N_forward(pp,sp,sf,cp,maxage,soil_mass,x,X); % mulitply concentration with (-1) to turn Maximum in nuclide concentratoin into minimum  for fminsearch

% run optimization to find maximum nuclide concentration
[D_Nmax,Nmax] = fminsearchbnd(fun,x0,D(1),D(2),options);


end

