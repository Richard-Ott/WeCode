function [XMAP,DMAP,W] = pairedCRN_Optim(pars10,pars36,D,X)
% Calculates the the "real" landscape denudation rate from a paired nuclide
% measurement and soil or bedrock chemistry data.
% Input: 
%           - pars10: Cronus parameters for 10Be computation (see 
%           Cronus_prep10.m)
%           - pars36: Cronus parameters for 36Cl computation
%           - D: [Dmin, Dmax] defines search boundaries for denudation rate
%           in mm/ka
%           - X: compositional structure (see CosmoDataRead.m)
% Output:   
%           - XMAP: compositional structure containing also the predicted
%           composition of the non-input parts (soil or bedrock)
%           - MAP: denudation rate in mm/ka
%           - WMAP: weathering rate in mm/ka
%
% Richard Ott, 2021

global scaling_model
v2struct(pars10);
v2struct(pars36);
soil_mass = X.soil_mass;

tic
% ------------------------------------------------------------------- %
% set some constants
sp10.depthtotop = soil_mass;             % set depth to soil bedrock interface
sp36.depthtotop = soil_mass;             % set depth to soil bedrock interface

% Figure out the maximum possible depth at which we'll ever need a
% production rate.  
Nobs = [nominal10(9);nominal36(1)];       % measured concentration

% PRIORS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% denudation
D = D/10*sp10.rb;                         % convert to g/cm²/ka for Cronus

% Qz concentration
switch X.mode
    case 'soil'
        fQz = [0,X.fQzS];                 % fraction of quartz in bedrock (cannot be smaller than fraction in soil)
    case 'bedrock'
        fQz = [X.fQzB,1-X.fXB];           % fraction of quartz in soil (cannot be smaller than fraction in bedrock)
end

% Resolution, under which parameters resolution does stop the model
tolerance = 0.005;                        % as fraction, I set it to 0.01 so that the total error of N10 and N36 is together never more than 1% of the nuclide concentration

% INIITAL MODEL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x0(1) = fQz(1)+rand*diff(fQz);            % create random quartz concentration

% For the dendation rate, the parameter seach starts at what would
% be the conventional 36Cl-denudation rate as a best guess.
spini = sp36; spini.depthtotop = 0;            % set depth to top = 0 for initial erosion rate guess
x0(2) = cl36erateraw(pp,spini,sf36,cp36,scaling_model,0);  % start parameters search at 36Cl rate (for a carbonate composition this should be closer to the real denudation rate)


%% OPTIMIZATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
options = optimset('MaxIter',5e4,'TolFun',tolerance,'TolX',0.05);            
% These options may need to be tuned specifically to your problem
% TolFun, maximum value that the function is allowed to be off, at/g
% TolX, the tolerance value in x-direction (erosion rates in g/cm2/a), which the
% function will use for subsequent iterations as tolerance level for
% stopping the algorithm. 

% define optimization function
fun = @(x) paired_N_forward_OptimWrapper(pp,sp10,sp36,sf10,sf36,cp10,cp36,maxage10,maxage36,scaling_model,soil_mass,x,X,Nobs);
[MAP,~] = fminsearchbnd(fun,x0,[fQz(1),D(1)],[fQz(2),D(2)],options);

%% get composition 
[~,~,DMAP, W, XMAP] = paired_N_forward(pp,sp10,sp36,sf10,sf36,cp10,cp36,maxage10,maxage36,scaling_model,soil_mass,MAP,X);

% convert denudation/weathering rates from g/cm2/ka to mm/ka
DMAP = DMAP/sp10.rb*10;  
W = W/sp10.rb*10;

end

