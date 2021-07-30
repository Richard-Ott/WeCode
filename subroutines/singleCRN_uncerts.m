function [MAP_uncerts, X_uncerts] = singleCRN_uncerts(pars,scaling,D,X,MAP,thres)
% Calculates the "real" denudation rate from a nuclide measurement, the
% bedrock or soil chemistry and a weathering rate.
% The solution is found through a MCMC algorithm.
% Richard Ott, 2021

wb = waitbar(0,'calculating uncertainties...');

%% NUCLIDE CONCENTRATION AND WEATHERING UNCERTAINTY
% only uncertainties from nuclide concentration and weathering rate are
% considered
uncertainty=0.0;
X_uncert = [0,0,0];
for i=1:2
    parsC = pars;
    if i == 1
        if X.n == 1
            thisdelta=0.1*pars.nominal10(9);
            parsC.nominal10(9) = pars.nominal10(9) +  thisdelta;
            parsC.sp10.concentration10 = parsC.nominal10(9);
            parsC.cp10.N36m = parsC.nominal10(9);
        elseif X.n == 2
            thisdelta=0.1*pars.nominal36(1);
            parsC.nominal36(1) = pars.nominal36(1) +  thisdelta;
            parsC.sp36.concentration36 = parsC.nominal36(1);
            parsC.cp36.N36m = parsC.nominal36(1);
        end
    else
      thisdelta=0.1*X.W;
    end
    
    % run inversion, try to start close to final value to speed up
    % inversion
    try
        [dX,deltaerate,~] = singleCRN_MCMC(parsC,scaling,[MAP- MAP/10, ;MAP+MAP/10],X,thres);
    catch
        [dX,deltaerate,~] = singleCRN_MCMC(parsC,scaling,D,X,thres);
    end
    
    deratei=(deltaerate - MAP)/(thisdelta);
    
    switch X.mode
        case 'soil'
            dXi = [dX.fQzB - X.fQzB, dX.fCaB - X.fCaB, dX.fXB - X.fXB] ./thisdelta;
        case 'bedrock'
            dXi = [dX.fQzS - X.fQzS, dX.fCaS - X.fCaS, dX.fXS - X.fXS] ./thisdelta;
    end
    
    if i == 1 % N uncertainty
        if X.n == 1
            uncertainty=uncertainty+(deratei^2*pars.uncerts10(9)^2);
            X_uncert    = dXi.^2 .* pars.uncerts10(9)^2;  % relative uncertainty in compositions
        elseif X.n == 2
            uncertainty=uncertainty+(deratei^2*pars.uncerts36(1)^2);
            X_uncert    = dXi.^2 .* pars.uncerts36(1)^2;            
        end 
        waitbar(0.33,wb)
    else     % W uncertainty
        uncertainty=uncertainty+(deratei^2*X.Wstd^2);
        X_uncert    = X_uncert + (dXi.^2 .* X.Wstd.^2);  % relative uncertainty in compositions
        waitbar(0.66,wb)
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PRODUCTION UNCERTAINTY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

deltapp = pars.pp;
if X.n == 1
    deltapp.PsBe= pars.pp.PsBe+0.1*abs(pars.pp.PsBe);
    % run inversion, try to start close to final value to speed up
    % inversion
    try
        [dX,deltaerate,~] = singleCRN_MCMC(pars,scaling,[MAP- MAP/10, ;MAP+MAP/10],X,thres);
    catch
        [dX,deltaerate,~] = singleCRN_MCMC(pars,scaling,D,X,thres);
    end
    deratepsBe  =(deltaerate-MAP)/(0.1*abs(pars.pp.PsBe));
    % compositional uncert
    switch X.mode
        case 'soil'
            dXi = [dX.fQzB - X.fQzB, dX.fCaB - X.fCaB, dX.fXB - X.fXB] ./thisdelta;
        case 'bedrock'
            dXi = [dX.fQzS - X.fQzS, dX.fCaS - X.fCaS, dX.fXS - X.fXS] ./thisdelta;
    end
    X_uncert    = X_uncert + (dXi.^2 .* pars.pp.sigmaPsBe.^2);  
    uncertainty = uncertainty+(deratepsBe^2*pars.pp.sigmaPsBe^2);
    
elseif X.n == 2
    deltapp.PsCa= pars.pp.PsCa0+0.1*abs(pars.pp.PsCa0);
    try
        [dX,deltaerate,~] = singleCRN_MCMC(pars,scaling,[MAP- MAP/10, ;MAP+MAP/10],X,thres);
    catch
        [dX,deltaerate,~] = singleCRN_MCMC(pars,scaling,D,X,thres);
    end
    
    % compositional uncert
    switch X.mode
        case 'soil'
            dXi = [dX.fQzB - X.fQzB, dX.fCaB - X.fCaB, dX.fXB - X.fXB] ./thisdelta;
        case 'bedrock'
            dXi = [dX.fQzS - X.fQzS, dX.fCaS - X.fCaS, dX.fXS - X.fXS] ./thisdelta;
    end
    
    X_uncert  = X_uncert + (dXi.^2 .* pars.pp.sigmaPsCa0.^2);
    deratepsCa  =(deltaerate-MAP)/(0.1*abs(pars.pp.PsCa0));
    uncertainty = uncertainty+(deratepsCa^2*pars.pp.sigmaPsCa0^2);    
end
MAP_uncerts     = sqrt(uncertainty);

switch X.mode
    case 'soil'
        X_uncerts       = sqrt(X_uncert).* [X.fQzB, X.fCaB, XfXB];
    case 'bedrock'
        X_uncerts       = sqrt(X_uncert).* [X.fQzS, X.fCaS, X.fXS];
end

close(wb)

end

