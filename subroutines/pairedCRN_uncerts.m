function [MAP_uncerts, X_uncerts] = pairedCRN_uncerts(pars10,pars36,D,X,XMAP,MAP)
% Richard Ott, 2021

wb = waitbar(0,'calculating uncertainties...');
global scaling_model

%% NUCLIDE CONCENTRATION UNCERTAINTY

parsC10 = pars10; parsC36 = pars36;
uncertainty=0.0;
X_uncert = [0,0,0];
for i=1:2
    if i == 1
        thisdelta=0.1*pars10.nominal10(9);                        % 10% of total N
        parsC10.nominal10(9) = pars10.nominal10(9) +  thisdelta;    % add 10% to total N
        parsC10.sp10.concentration10 = parsC10.nominal10(9);        % update sample parameters
        parsC10.cp10.N10m = parsC10.nominal10(9);                   % update computed parameters 
        parsC10.uncertFlag =1;
    else
        thisdelta=0.1*pars36.nominal36(1);
        parsC36.nominal36(1) = pars36.nominal36(1) +  thisdelta;
        parsC36.sp36.concentration36 = parsC36.nominal36(1);
        parsC36.cp36.N36m = parsC36.nominal36(1);
        parsC36.uncertFlag =1;
    end
    
    % run inversion, try to start close to final value to speed up
    % inversion, if it does not work-use full range
    try
        [dX,deltaerate,~] = pairedCRN_MCMC(parsC10,parscC36,[(2)- MAP(2)/10, ;MAP(2)+MAP(2)/10],X); % find erate and composition
    catch
        [dX,deltaerate,~] = pairedCRN_MCMC(parsC10,parsC36,D,X);
    end
    
    deratei=(deltaerate(2) - MAP(2))/thisdelta;
    
    switch X.mode
        case 'soil'
            dXi = [dX.fQzB - XMAP.fQzB, dX.fCaB - XMAP.fCaB, dX.fXB - XMAP.fXB] ./thisdelta;
        case 'bedrock'
            dXi = [dX.fQzS - XMAP.fQzS, dX.fCaS - XMAP.fCaS, dX.fXS - XMAP.fXS] ./thisdelta;
    end
    
    if i == 1 % N10 uncertainty
        uncertainty=uncertainty+ deratei^2*pars10.uncerts10(9)^2;
        X_uncert    = X_uncert + dXi.^2 .* pars10.uncerts10(9)^2;  % relative uncertainty in compositions
    else     % N36 uncertainty
        uncertainty=uncertainty+ deratei^2*pars36.uncerts36(1)^2;
        X_uncert    = X_uncert + dXi.^2 .* pars36.uncerts36(1)^2;    
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
        parsC10.pp = deltapp; parsC10.uncertFlag =1;
    else
        deltapp = pars36.pp;
        thisdelta = 0.1*abs(pars36.pp.PsCa0);
        deltapp.PsCa0= pars36.pp.PsCa0+ thisdelta;
        deltacp=comppars36(deltapp,pars36.sp36,pars36.sf36,pars36.maxdepth36);% new cp
        parsC36 = pars36;
        parsC36.cp36 = deltacp;
        parsC36.pp = deltapp; parsC36.uncertFlag =1;
    end
    
    % run inversion, try to start close to final value to speed up
    % inversion, if it does not work-use full range
    try
        [dX,deltaerate,~] = pairedCRN_MCMC(parsC10,parsC36,[MAP(2)- MAP(2)/10, ;MAP(2)+MAP(2)/10],X);
    catch
        [dX,deltaerate,~] = pairedCRN_MCMC(parsC10,parsC36,D,X);
    end
    
    derateps  =(deltaerate(2)-MAP(2))/ thisdelta;
    % compositional uncert
    switch X.mode
        case 'soil'
            dXi = [dX.fQzB - XMAP.fQzB, dX.fCaB - XMAP.fCaB, dX.fXB - XMAP.fXB] ./thisdelta;
        case 'bedrock'
            dXi = [dX.fQzS - XMAP.fQzS, dX.fCaS - XMAP.fCaS, dX.fXS - XMAP.fXS] ./thisdelta;
    end
    
    % add uncertainties up
    if i == 1
        X_uncert    = X_uncert +   dXi.^2 .* pars.pp.sigmaPsBe.^2;  
        uncertainty = uncertainty+ derateps^2*pars10.pp.sigmaPsBe^2;
    else
        X_uncert    = X_uncert +   dXi.^2 .* pars.pp.sigmaPsCa0.^2;  
        uncertainty = uncertainty+(deratepsCa^2*pars36.pp.sigmaPsCa0^2);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

