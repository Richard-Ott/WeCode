function D_uncerts = SBW_uncerts(pars,X,D)
% Richard Ott, 2021

wb = waitbar(0,'calculating uncertainties...');
global scaling_model


%% NUCLIDE CONCENTRATION AND WEATHERING UNCERTAINTY
% only uncertainties from nuclide concentration and weathering rate are
% considered in this section
uncertainty=0.0;
for i=1:2
    Xcur = X;
    parsC = pars;
    if i == 1  % first calculate uncertainty from N
        if Xcur.n == 1
            thisdelta=0.05*pars.nominal10(9);                       % 5% of total N
            parsC.nominal10(9) = pars.nominal10(9) +  thisdelta;    % add 5% to total N
            parsC.sp10.concentration10 = parsC.nominal10(9);        % update sample parameters
            parsC.cp10.N10m = parsC.nominal10(9);                   % update computed parameters 
        elseif Xcur.n == 2
            thisdelta=0.05*pars.nominal36(1);                  % for the soluble mineral and weathering I rather lower the value because otherwise you quickly run into a no-go zone
            parsC.nominal36(1) = pars.nominal36(1) +  thisdelta;
            parsC.sp36.concentration36 = parsC.nominal36(1);
            parsC.cp36.N36m = parsC.nominal36(1);
        end
    else
      thisdelta=0.05*Xcur.W;
      Xcur.W = Xcur.W - thisdelta;
    end
    
    % run optimization
    deltaerate = SBW_search(parsC,Xcur); % find erate and composition

    deratei=(deltaerate - D)/thisdelta;
    
    if i == 1 % N uncertainty
        if X.n == 1
            uncertainty=uncertainty+ deratei^2*pars.uncerts10(9)^2;
        elseif X.n == 2
            uncertainty=uncertainty+ deratei^2*pars.uncerts36(1)^2;
        end 
        waitbar(0.33,wb)
    else     % W uncertainty
        uncertainty=uncertainty+ deratei^2*X.Wstd^2;
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
    parsC.pp = deltapp; 
    
    % run optimization
    deltaerate = SBW_search(parsC,X);

    deratepsBe  =(deltaerate-D)/thisdelta;

    uncertainty = uncertainty+ deratepsBe^2*pars.pp.sigmaPsBe^2;
    
elseif X.n == 2
    thisdelta = 0.1*abs(pars.pp.PsCa0);
    deltapp.PsCa= pars.pp.PsCa0+thisdelta;
    
    [deltaerate,dX] = SBW_search(parsC,X);
    deratepsCa  =(deltaerate-D)/thisdelta;
    uncertainty = uncertainty+(deratepsCa^2*pars.pp.sigmaPsCa0^2);  
    waitbar(1,wb)
end

%Combined uncertainty of denudation rate
D_uncerts     = sqrt(uncertainty);

close(wb)

end

