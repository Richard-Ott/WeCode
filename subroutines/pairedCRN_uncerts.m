function [MAP_uncerts, X_uncerts, W_uncert] = pairedCRN_uncerts(pars10,pars36,D,X,XMAP,MAP,WMAP)
% Calculates the uncerttainty in denudation rate, composition, and
% weathering rate from paired nuclide measurements.
% This script only includes uncertainties that arise from the two nuclide
% concentrations, and their spallation production factors.
% Input: 
%       - pars10 & pars 36: Cronus parameters for CRN calculation see
%       (Cronus_prepXX.m)
%       - D: [Dmin, Dmax] defines search boundaries for denudation rate
%           in mm/ka
%       - X: compositional structure (see CosmoDataRead.m)
%       - XMAP: best fit compositon
%       - MAP: best-fit denudation rate (mm/ka)
%       - WMAP: best-fit weathering rate (mm/ka)
% Output: 
%       - MAP_uncerts: Uncertainty denudation rate (mm/ka)
%       - X_uncerts:   Uncertainty composition 
%       - W_uncerts:   Uncertainty weathering rate (mm/ka)
%
% Richard ott, 2021

wb = waitbar(0,'calculating uncertainties...');
global scaling_model

%% NUCLIDE CONCENTRATION UNCERTAINTY

uncertainty=0;
X_uncert   = [0,0,0];
W_uncert   = 0;
for i=1:2
    parsC10 = pars10; parsC36 = pars36;
    if i == 1
        thisdelta=0.05*pars10.nominal10(9);                          % 10% of total N
        parsC10.nominal10(9) = pars10.nominal10(9) +  thisdelta;    % add 10% to total N
        parsC10.sp10.concentration10 = parsC10.nominal10(9);        % update sample parameters
        parsC10.cp10.N10m = parsC10.nominal10(9);                   % update computed parameters 
    else
        thisdelta=0.05*pars36.nominal36(1);
        parsC36.nominal36(1) = pars36.nominal36(1) +  thisdelta; % if W is close to D, then the addition might lead to problems and one should acutally subtract here
        parsC36.sp36.concentration36 = parsC36.nominal36(1);
        parsC36.cp36.N36m = parsC36.nominal36(1);
    end
    
    [dX,deltaerate,deltaW] = pairedCRN_Optim(parsC10,parsC36,D,X);

    % find differences
    deratei =(deltaerate - MAP)/thisdelta;
    deltaWi =(deltaW - WMAP)/thisdelta;
    switch X.mode
        case 'soil'
            dXi = [dX.fQzB - XMAP.fQzB, dX.fCaB - XMAP.fCaB, dX.fXB - XMAP.fXB] ./thisdelta;
        case 'bedrock'
            dXi = [dX.fQzS - XMAP.fQzS, dX.fCaS - XMAP.fCaS, dX.fXS - XMAP.fXS] ./thisdelta;
    end
    
    if i == 1 % N10 uncertainty
        uncertainty=uncertainty+ deratei^2*pars10.uncerts10(9)^2;
        X_uncert    = X_uncert + dXi.^2  .*pars10.uncerts10(9)^2;  % relative uncertainty in compositions
        W_uncert    = W_uncert + deltaWi^2*pars10.uncerts10(9)^2;
        waitbar(0.25,wb)
    else     % N36 uncertainty
        uncertainty=uncertainty+ deratei^2*pars36.uncerts36(1)^2;
        X_uncert    = X_uncert + dXi.^2 .* pars36.uncerts36(1)^2;    
        W_uncert    = W_uncert + deltaWi^2*pars36.uncerts36(1)^2;
        waitbar(0.5,wb)
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PRODUCTION UNCERTAINTY (only spallation) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:2
    if i == 1
        deltapp = pars10.pp;
        thisdelta = 0.1*abs(pars10.pp.PsBe);
        deltapp.PsBe= pars10.pp.PsBe+ thisdelta;  % 10% change in spallation production
        deltacp=comppars1026(deltapp,pars10.sp10,pars10.sf10,pars10.maxdepth10);% new cp
        parsC10 = pars10;
        parsC10.cp10 = deltacp;
        parsC10.pp = deltapp; 
    else
        deltapp = pars36.pp;
        thisdelta = 0.1*abs(pars36.pp.PsCa0);
        deltapp.PsCa0= pars36.pp.PsCa0+ thisdelta;
        deltacp=comppars36(deltapp,pars36.sp36,pars36.sf36,pars36.maxdepth36);% new cp
        parsC36 = pars36;
        parsC36.cp36 = deltacp;
        parsC36.pp = deltapp; 
    end
    

    [dX,deltaerate,deltaW] = pairedCRN_Optim(parsC10,parsC36,D,X);
    
    % differences
    derateps  =(deltaerate-MAP)/ thisdelta;
    deltaWi   = (deltaW - WMAP)/thisdelta;
    % compositional uncert
    switch X.mode
        case 'soil'
            dXi = [dX.fQzB - XMAP.fQzB, dX.fCaB - XMAP.fCaB, dX.fXB - XMAP.fXB] ./thisdelta;
        case 'bedrock'
            dXi = [dX.fQzS - XMAP.fQzS, dX.fCaS - XMAP.fCaS, dX.fXS - XMAP.fXS] ./thisdelta;
    end
    
    % add uncertainties up
    if i == 1
        X_uncert    = X_uncert +   dXi.^2   .* pars10.pp.sigmaPsBe^2;  
        uncertainty = uncertainty+ derateps^2* pars10.pp.sigmaPsBe^2;
        W_uncert    = W_uncert +   deltaWi^2 * pars10.pp.sigmaPsBe^2;
        waitbar(0.75,wb)
    else
        X_uncert    = X_uncert +   dXi.^2  .* pars36.pp.sigmaPsCa0^2;  
        uncertainty = uncertainty+derateps^2* pars36.pp.sigmaPsCa0^2;
        W_uncert    = W_uncert +   deltaWi^2* pars36.pp.sigmaPsCa0^2;
        waitbar(1,wb)
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Combined uncertainty of denudation rate
MAP_uncerts     = sqrt(uncertainty);
W_uncert        = sqrt(W_uncert);

% Combined uncertainty of mineral fractions
switch X.mode
    case 'soil'
        X_uncerts       = sqrt(X_uncert);
    case 'bedrock'
        X_uncerts       = sqrt(X_uncert);
end

close(wb)

end

