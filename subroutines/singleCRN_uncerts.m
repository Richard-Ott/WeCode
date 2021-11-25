function MAP_uncerts = singleCRN_uncerts(pars,D,X,MAP,thres)
% Richard Ott, 2021

wb = waitbar(0,'calculating uncertainties...');
global scaling_model

bounds = D;   


%% NUCLIDE CONCENTRATION AND CDF UNCERTAINTY
% only uncertainties from nuclide concentration and weathering rate are
% considered in this section
uncertainty=0.0;
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
            if pars.nonunique && pars.nominal36(1)+thisdelta > pars.Nmax
                parsC.nominal36(1) = pars.nominal36(1) -  thisdelta;
            else
                parsC.nominal36(1) = pars.nominal36(1) +  thisdelta;
            end
            parsC.sp36.concentration36 = parsC.nominal36(1);
            parsC.cp36.N36m = parsC.nominal36(1);
        end
    else
      thisdelta=0.05*Xcur.CDF;
      Xcur.CDF = Xcur.CDF - thisdelta;
    end
    
    % run optimization
    deltaerate = singleCRN_Optim(parsC,bounds,Xcur,thres); % find erate and composition
    if isfield(pars,'nonunique'); [deltaerate,ind] = max(deltaerate); dX = dX(ind);end

    deratei=(deltaerate - MAP)/thisdelta;
    
    if i == 1 % N uncertainty
        if X.n == 1
            uncertainty=uncertainty+ deratei^2*pars.uncerts10(9)^2;
        elseif X.n == 2
            uncertainty=uncertainty+ deratei^2*pars.uncerts36(1)^2;
        end 
        waitbar(0.33,wb)
    else     % CDF uncertainty
        uncertainty=uncertainty+ deratei^2*X.CDFstd^2;
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
    
    % run optimization
    deltaerate = singleCRN_Optim(parsC,bounds,X,thres);

    deratepsBe  =(deltaerate-MAP)/thisdelta;

    uncertainty = uncertainty+ deratepsBe^2*pars.pp.sigmaPsBe^2;
    
elseif X.n == 2
    thisdelta = 0.1*abs(pars.pp.PsCa0);
    deltapp.PsCa= pars.pp.PsCa0+thisdelta;
    
    
   deltaerate = singleCRN_Optim(pars,bounds,X,thres);
    if isfield(pars,'nonunique'); [deltaerate,ind] = max(deltaerate); dX = dX(ind);end

    deratepsCa  =(deltaerate-MAP)/thisdelta;
    uncertainty = uncertainty+(deratepsCa^2*pars.pp.sigmaPsCa0^2);  
    waitbar(1,wb)
end

%Combined uncertainty of denudation rate
MAP_uncerts     = sqrt(uncertainty);
close(wb)

end

