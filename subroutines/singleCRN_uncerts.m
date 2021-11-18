function [MAP_uncerts, X_uncerts] = singleCRN_uncerts(pars,D,X,XMAP,MAP,thres)
% Richard Ott, 2021

wb = waitbar(0,'calculating uncertainties...');
global scaling_model

%% NUCLIDE CONCENTRATION AND WEATHERING UNCERTAINTY
% only uncertainties from nuclide concentration and weathering rate are
% considered
uncertainty=0.0;
X_uncert = [0,0,0];
for i=1:2
    parsC = pars;
    parsC.uncertFlag =1;
    if i == 1
        if X.n == 1
            thisdelta=0.1*pars.nominal10(9);                       % 1% of total N
            parsC.nominal10(9) = pars.nominal10(9) +  thisdelta;    % add 1% to total N
            parsC.sp10.concentration10 = parsC.nominal10(9);        % update sample parameters
            parsC.cp10.N10m = parsC.nominal10(9);                   % update computed parameters 
        elseif X.n == 2
            thisdelta=0.1*pars.nominal36(1);                        % for the soluble mineral and weathering I rather lower the value because otherwise you quickly run into a no-go zone
            parsC.nominal36(1) = pars.nominal36(1) -  thisdelta;
            parsC.sp36.concentration36 = parsC.nominal36(1);
            parsC.cp36.N36m = parsC.nominal36(1);
        end
    else
      thisdelta=0.1*X.W;
      X.W = X.W - thisdelta;
    end
    
    % run inversion, try to start close to final value to speed up
    % inversion, if it does not work-use full range
    try
        [dX,deltaerate,~] = singleCRN_MCMC(parsC,[MAP- MAP/2, ;MAP+MAP/2],X,thres); % find erate and composition
    catch
        [dX,deltaerate,~] = singleCRN_MCMC(parsC,D,X,thres);
    end
    
    deratei=(deltaerate - MAP)/thisdelta;
    
    switch X.mode
        case 'soil'
            dXi = [dX.fQzB - XMAP.fQzB, dX.fCaB - XMAP.fCaB, dX.fXB - XMAP.fXB] ./thisdelta;
        case 'bedrock'
            dXi = [dX.fQzS - XMAP.fQzS, dX.fCaS - XMAP.fCaS, dX.fXS - XMAP.fXS] ./thisdelta;
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
    deltapp.PsBe= pars.pp.PsBe+ thisdelta;  % 1% chnage in spallation production
    % run inversion, try to start close to final value to speed up
    % inversion
    parsC = pars;
    parsC.pp = deltapp; parsC.uncertFlag =1;
    % run inversion, try to start close to final value to speed up
    % inversion, if it does not work-use full range
    try
        [dX,deltaerate,~] = singleCRN_MCMC(parsC,[MAP- MAP/2, ;MAP+MAP/2],X,thres);
    catch
        [dX,deltaerate,~] = singleCRN_MCMC(parsC,D,X,thres);
    end
    deratepsBe  =(deltaerate-MAP)/thisdelta;
    % compositional uncert
    switch X.mode
        case 'soil'
            dXi = [dX.fQzB - XMAP.fQzB, dX.fCaB - XMAP.fCaB, dX.fXB - XMAP.fXB] ./thisdelta;
        case 'bedrock'
            dXi = [dX.fQzS - XMAP.fQzS, dX.fCaS - XMAP.fCaS, dX.fXS - XMAP.fXS] ./thisdelta;
    end
    X_uncert    = X_uncert +   dXi.^2 .* pars.pp.sigmaPsBe.^2;  
    uncertainty = uncertainty+ deratepsBe^2*pars.pp.sigmaPsBe^2;
    
elseif X.n == 2
    thisdelta = 0.1*abs(pars.pp.PsCa0);
    deltapp.PsCa= pars.pp.PsCa0+thisdelta;
    try
        [dX,deltaerate,~] = singleCRN_MCMC(pars,[MAP- MAP/10, ;MAP+MAP/10],X,thres);
    catch
        [dX,deltaerate,~] = singleCRN_MCMC(pars,D,X,thres);
    end
    
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

close(wb)

end

