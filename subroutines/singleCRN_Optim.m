function varargout = singleCRN_Optim(pars,D,X,thres)
% Calculates the "real" denudation rate from a nuclide measurement, the
% bedrock or soil chemistry, and a weathering rate.
% The solution is found through the fminsearch optimization algorithm.
% % Input: 
%           - pars: Cronus parameters for CRN computation (see 
%           Cronus_prepXX.m)
%           - D: [Dmin, Dmax] defines search boundaries for denudation rate
%           in mm/ka
%           - X: compositional structure (see CosmoDataRead.m)
%           - thres: threshold of nuclide concentration error at which 
%           optimization has converged in % of N
% Output:   
%           - X_MAP: compositional structure containing also the predicted
%           composition of the non-input parts (soil or bedrock)
%           - MAP: denudation rate in mm/ka
%           --> only MAP output for X.mode = noComp (no composition defined)
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
sp.depthtotop = soil_mass;           % set depth to soil bedrock interface
% Figure out the maximum possible depth at which we'll ever need a
% production rate.  
data_ind = n+8-9*(n-1);                % index I use for referencing into the nominal and uncerts arrays, because unfortunately Cronus uses different order for the data depending on the nuclide
Nobs = nominal(data_ind);              % measured concentration
        
% INPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Minimum denudation rate %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The enrichment cannot produce denudation rates that would lead to
% enrichment of insoluble minerals that go beyonf 100% of the soil fraction
switch X.mode; case 'bedrock'; insol = X.fQzB + X.fXB; case 'noComp'; insol= -1; end
if exist('insol','var'); D(1) = W/(1-insol); end


D = D/10*sp.rb;                        % convert to g/cm²/ka for Cronus
W = W/10*sp.rb;                        % convert to g/cm²/ka 
tolerance = Nobs*thres/100;            % in at/g for nuclides

% Define options for optimization  
options = optimset('MaxIter',5e4,'TolFun',tolerance);            
% These options may need to be tuned specifically to your problem
% TolFun, maximum value that the function is allowed to be off, at/g
% TolX, the tolerance value in x-direction (erosion rates in g/cm2/a), which the
% function will use for subsequent iterations as tolerance level for
% stopping the algorithm. 

x0model = {@be10erateraw, @cl36erateraw};
spini = sp; spini.depthtotop = 0;            % set depth to top = 0 for initial erosion rate guess
x0 = x0model{n}(pp,spini,sf,cp,scaling_model,0);
if strcmpi(X.mode,'noComp')
    fun = @(x) abs(N10_forward_noComp(pp,sp,sf,cp,maxage,scaling_model,soil_mass,x,X) - Nobs);
else
    fun = @(x) abs(Comp_and_N_forward(pp,sp,sf,cp,maxage,soil_mass,x,X) - Nobs);
end
[MAP,~] = fminsearchbnd(fun,x0,D(1),D(2),options);

if MAP == D(1)
    warning('The best-fit denudation rate equals the minimum denudation rate. Check your input and think about your assumptions.')
    [nobound_MAP,~] = fminsearchbnd(fun,x0,0,D(2),options);
    disp(['Without input constraints the best-fit denudation rate would be ' num2str(nobound_MAP/sp.rb*10) ' mm/ka'])
end

toc

% for soluble mineral with a nuclide concentration that is nonunique, we
% need to calculate the second solution.
if exist('nonunique')
    if MAP > D_Nmax; bounds = [D(1), D_Nmax]; else; bounds = [D_Nmax,D(2)];end
    fun = @(x) abs(Comp_and_N_forward(pp,sp,sf,cp,maxage,soil_mass,x,X) - Nobs);
    MAP = [MAP, fminsearchbnd(fun,x0,bounds(1),bounds(2),options)];
end    

MAP = MAP/sp.rb*10;  % convert back to mm/ka

varargout{1} = MAP;
if nargout == 2
    % run forward model one more time for compositional output
    [~,X_MAP] = Comp_and_N_forwardX(pp,sp,sf,cp,maxage,soil_mass,MAP(1)/10*sp.rb,X);
    if exist('nonunique')
        [~,X_MAP(2)] = Comp_and_N_forwardX(pp,sp,sf,cp,maxage,soil_mass,MAP(2)/10*sp.rb,X);
    end

    varargout{2} = X_MAP;
end
end

