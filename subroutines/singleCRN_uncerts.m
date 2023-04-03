function varargout = singleCRN_uncerts(pars,D,X,MAP,thres,varargin)
% Calculates the uncertainty in denudation rate and composition from a
% nuclide measurement with weathering rate.
% This script only includes uncertainties that arise from the nuclide
% concentrations, weathering rate, and the spallation production factor.
% Input: 
%       - pars: Cronus parameters for CRN calculation see (Cronus_prepXX.m)
%       - D: [Dmin, Dmax] defines search boundaries for denudation rate
%           in mm/ka
%       - X: compositional structure (see CosmoDataRead.m
%       - MAP: best-fit denudation rate (mm/ka)
%       - WMAP: best-fit weathering rate (mm/ka)
%       - OPTIONAL: XMAP -> best fit compositon
% Output: 
%       - MAP_uncerts: Uncertainty denudation rate (mm/ka)
%       - X_uncerts:   Uncertainty composition 
%       --> only MAP_uncerts for X.mode == noComp (no composition defined)
% Richard ott, 2021

if nargin == 6
    XMAP = varargin{1};
end
wb = waitbar(0,'calculating uncertainties...');
global scaling_model

% if there were two solutions, only calculate uncertainty for the higher
% denudation rate solution
if isfield(pars,'nonunique')
    [MAP,ind] = max(MAP); XMAP = XMAP(ind);
    bounds = [pars.D_Nmax, ;D(2)];  % bounds for upper rate
else
    bounds = D;  
end    


%% NUCLIDE CONCENTRATION AND WEATHERING UNCERTAINTY
% only uncertainties from nuclide concentration and weathering rate are
% considered in this section
uncertainty=0.0;
X_uncert = [0,0,0];
for i=1:2
    Xcur = X;
    parsC = pars;
    parsC.uncertFlag =1;
    if i == 1  % first calculate uncertainty from N
        if Xcur.n == 1
            thisdelta=0.05*pars.nominal10(9);                       % 5% of total N
            parsC.nominal10(9) = pars.nominal10(9) +  thisdelta;    % add 5% to total N
            parsC.sp10.concentration10 = parsC.nominal10(9);        % update sample parameters
            parsC.cp10.N10m = parsC.nominal10(9);                   % update computed parameters 
        elseif Xcur.n == 2
            thisdelta=0.05*pars.nominal36(1);                  % for the soluble mineral and weathering I rather lower the value because otherwise you quickly run into a no-go zone
            if isfield(pars,'nonunique') && pars.nominal36(1)+thisdelta > pars.Nmax
                parsC.nominal36(1) = pars.nominal36(1) -  thisdelta;
            else
                parsC.nominal36(1) = pars.nominal36(1) +  thisdelta;
            end
            parsC.sp36.concentration36 = parsC.nominal36(1);
            parsC.cp36.N36m = parsC.nominal36(1);
        end
    else
      thisdelta=0.05*Xcur.W;
      Xcur.W = Xcur.W - thisdelta;
    end
    
    % run optimization
    if strcmpi(X.mode,'noComp')
        deltaerate = singleCRN_Optim(parsC,bounds,Xcur,thres); % find erate and composition
    else
        [deltaerate,dX] = singleCRN_Optim(parsC,bounds,Xcur,thres); % find erate and composition
    end
    if isfield(pars,'nonunique'); [deltaerate,ind] = max(deltaerate); dX = dX(ind);end

    deratei=(deltaerate - MAP)/thisdelta;
    
    switch X.mode
        case 'soil'
            dXi = [dX.fQzB - XMAP.fQzB, dX.fCaB - XMAP.fCaB, dX.fXB - XMAP.fXB] ./thisdelta;
        case 'bedrock'
            dXi = [dX.fQzS - XMAP.fQzS, dX.fCaS - XMAP.fCaS, dX.fXS - XMAP.fXS] ./thisdelta;
        case 'noComp'
            dXi = [0,0,0];  % dummy variables
    end
    
    if i == 1 % N uncertainty
        if X.n == 1
            uncertainty=uncertainty+ deratei^2*pars.uncerts10(9)^2;
            X_uncert    = X_uncert + dXi.^2 .* pars.uncerts10(9)^2;  % relative uncertainty in compositions
        elseif X.n == 2
            uncertainty=uncertainty+ deratei^2*pars.uncerts36(1)^2;
            X_uncert    = X_uncert + dXi.^2 .* pars.uncerts36(1)^2;            
        end 
        waitbar(0.33,wb)
    else     % W uncertainty
        uncertainty=uncertainty+ deratei^2*X.Wstd^2;
        X_uncert    = X_uncert + dXi.^2 .* X.Wstd.^2;  % relative uncertainty in compositions
        waitbar(0.66,wb)
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PRODUCTION UNCERTAINTY (only spallation) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

deltapp = pars.pp;
if X.n == 1
    thisdelta = 0.1*abs(pars.pp.PsBe);
    deltapp.PsBe= pars.pp.PsBe+ thisdelta;  % 10% chnage in spallation production
    % run inversion, try to start close to final value to speed up
    % inversion
    parsC = pars;
    parsC.pp = deltapp; parsC.uncertFlag =1;
    
    % run optimization
    if strcmpi(X.mode,'noComp')
        deltaerate = singleCRN_Optim(parsC,bounds,X,thres);
    else
        [deltaerate,dX] = singleCRN_Optim(parsC,bounds,X,thres);
    end

    deratepsBe  =(deltaerate-MAP)/thisdelta;
    % compositional uncert
    switch X.mode
        case 'soil'
            dXi = [dX.fQzB - XMAP.fQzB, dX.fCaB - XMAP.fCaB, dX.fXB - XMAP.fXB] ./thisdelta;
        case 'bedrock'
            dXi = [dX.fQzS - XMAP.fQzS, dX.fCaS - XMAP.fCaS, dX.fXS - XMAP.fXS] ./thisdelta;
        case 'noComp'
            dXi = [0,0,0];  % dummy variables
    end
    X_uncert    = X_uncert +   dXi.^2 .* pars.pp.sigmaPsBe.^2;  
    uncertainty = uncertainty+ deratepsBe^2*pars.pp.sigmaPsBe^2;
    
elseif X.n == 2
    thisdelta = 0.1*abs(pars.pp.PsCa0);
    deltapp.PsCa= pars.pp.PsCa0+thisdelta;
    
    
   [deltaerate,dX] = singleCRN_Optim(pars,bounds,X,thres);
    if isfield(pars,'nonunique'); [deltaerate,ind] = max(deltaerate); dX = dX(ind);end
    
    % compositional uncert
    switch X.mode
        case 'soil'
            dXi = [dX.fQzB - XMAP.fQzB, dX.fCaB - XMAP.fCaB, dX.fXB - XMAP.fXB] ./thisdelta;
        case 'bedrock'
            dXi = [dX.fQzS - XMAP.fQzS, dX.fCaS - XMAP.fCaS, dX.fXS - XMAP.fXS] ./thisdelta;
    end
    
    X_uncert  = X_uncert + dXi.^2 .* pars.pp.sigmaPsCa0.^2;
    deratepsCa  =(deltaerate-MAP)/thisdelta;
    uncertainty = uncertainty+(deratepsCa^2*pars.pp.sigmaPsCa0^2);  
    waitbar(1,wb)
end

%Combined uncertainty of denudation rate
MAP_uncerts     = sqrt(uncertainty);

% Combined uncertainty of mineral fractions
switch X.mode
    case 'soil'
        X_uncerts       = sqrt(X_uncert);
    case 'bedrock'
        X_uncerts       = sqrt(X_uncert);
end

varargout{1} = MAP_uncerts;
if nargout == 2
    varargout{2} = X_uncerts;
end
    

close(wb)

end

