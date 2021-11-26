function pars = solCRN_check(pars,D,X,err,plotFlag)
% Checks if the nuclide concentration for the soluble mineral is below the
% concentration for Dmin, or below the maximum possible concentration Nmax.
% Input: 
%       - pars: Cronus parameters for CRN calculation see (Cronus_prepXX.m)
%       - D: [Dmin, Dmax] defines search boundaries for denudation rate
%           in mm/ka
%       - X: compositional structure (see CosmoDataRead.m)
%       - err: threshold of nuclide concentration error at which 
%           optimization has converged in % of N
%       - plotFlag: 1 or 0, do you want to see a plot of nuclide
%       concentration vs. denudation rate. This helps you to see how close
%       you are to the zone of 2 solutions for N36.
%
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

% compute nuclide concentration for Dmin. If observed concentration is
% lower than for Dmin, then the calculation will have a unique solution
[N_Dmin,~] = Comp_and_N_forwardX(pp,sp,sf,cp,maxage,soil_mass,D(1),X);

if Nobs < N_Dmin
    sprintf('Your nuclide concentration is lower as for Dmin. All good.\nThere is a unique denudation rate for your measurement')
else
    sprintf('Your nuclide concentration is higher as for Dmin.\nThere will be 2 solutions to your denudation rate.\nThe code will now find the minimum denudation rate for a unique solution.')
    pars.nonunique = true;
end

%% Compute nuclide concentrations for a range of denudation rates
% provide visual if desired
if plotFlag
    Drates = linspace(D(1),D(2));   % denudation rates to compute N for
    Nm = arrayfun(@(x) Comp_and_N_forward(pp,sp,sf,cp,maxage,soil_mass,x,X), Drates);

    plot(Drates./sp.rb*10,Nm,'k-')  % N-modelled
    hold on 
    yline(Nobs,'r--');              % measured nuclide concentration   
    ylabel('nuclide concentration at/g')
    xlabel('Denudation rate mm/ka')
end

%% Identify maximum possible nuclide concentration for this sample and model
tolerance = Nobs*err/100;            % in at/g for nuclides
% Define options for optimization  
options = optimset('MaxIter',5e4,'TolFun',tolerance,'TolX',0.05); 

x0 = D(1)+3;       % start at denudation rate slightly higher than minimum
fun = @(x) (-1)*Comp_and_N_forward(pp,sp,sf,cp,maxage,soil_mass,x,X); % mulitply concentration with (-1) to turn Maximum in nuclide concentratoin into minimum  for fminsearch

% run optimization to find maximum nuclide concentration
[D_Nmax,Nmax] = fminsearchbnd(fun,x0,D(1),D(2),options);

if Nobs > Nmax*(-1)
    error('Your input concentration is higher as can be explained by the model parameters. Please, check input')
else
    sprintf(['The maximum modelled nuclide concentration is found for a denudation rate of ' num2str(D_Nmax/sp.rb*10) ' mm/ka \nIf you are confident that the denudation of your sample is higher, then there is a unique solution.'])
    pars.D_Nmax = D_Nmax;
    pars.Nmax = Nmax*(-1);
end


end

