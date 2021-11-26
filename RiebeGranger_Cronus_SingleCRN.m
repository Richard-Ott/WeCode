% The script shows how to calculate the denudation rate from a single nuclide
% measurement of a soluble or an insoluble target mineral. 
%
% The parameter search is performed via a Markov-Chain Monte Carlo approach
% with a Metropolis Hastings sampling alrogithm. 
%
% Richard Ott, 2021
%
% the current version is written for 10Be and 36Cl, could easily be
% expanded to other nuclides

clc
clear 
close all
addpath '.\subroutines'

% load data
[num,sampName,X,DEMdata] = CosmoDataRead('Test_Input_Single_Be.xlsx');

%% assign data and initial basin calculations --------------------------- %

% in case your denudation rate is from an alluvial sample and you desire a
% pixel-by-pixel calculated production rate provide a DEM
% THIS STEP REQUIRES TOPOTOOLBOX FUNCTIONS (Schwanghart & Scherler, 2014)
if strcmpi('basin',DEMdata.method)
    DEMdata.DEM = GRIDobj();        % interactively choose the DEM that encompasses the basin
    DEMdata.export = 1;             % do you want to save the individual sample scaling factors as .mat file?
    % This can be useful when the computation for scaling schemes like 'sa'
    % and 'sf'  takes a long time for a big basin and you want the scaling
    % factors saved for later
    
    [DEMdata.DB,DEMdata.utmzone] = getBasins(DEMdata.DEM,num(:,2),num(:,1),'ll');  % delineate drainage basins and check their geometry
end


%% Calculate production rates ------------------------------------------- %

Cronus_prep = {@Cronus_prep10, @Cronus_prep36};

pars = Cronus_prep{X.n}(num,DEMdata);

%% Run Optimization for "real" denudation rate -------------------------- %
% PRIORS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
D = [50,1e3];                          % Denudation min/max in mm/ka (Dmin > Weathering rate)
thres = 0.1;                           % threshold of nuclide concentration error at which optimization has converged in % of N

% if you're looking at the soluble target mineral, there might be more than
% one solution. This function will check for you if the parameters you
% suppied might lead to a nonunique solution. It's not necessary to run
% this, but advisable.
if X.n == 2 
    pars = solCRN_check(pars,D,X,thres,0);
end

% run inversion
[MAP,X_MAP] = singleCRN_Optim(pars,D,X,thres);   

% OPTIONAL
% estimate uncertainty, this take quite some time to calculate, therefore I
% commented the next line to speed up the example
[MAP_uncerts, X_uncerts] = singleCRN_uncerts(pars,D,X,X_MAP,MAP,thres);

%% OUTPUT RESULTS ------------------------------------------------------- %

% outputs the values for the higher rate if there was a nonunique solution
try
    % disp(['Denudation rate = ' num2str(round(MAP(end))) ' mm/ka'])
    disp(['Denudation rate = ' num2str(round(MAP(end))) ' ' char(177) ' ' num2str(round(MAP_uncerts)) ' mm/ka'])
    switch X.mode
        case 'soil'
            disp(['Fraction of quartz in bedrock fQzB = ' num2str(X_MAP(end).fQzB) ' ' char(177) ' ' num2str(round(X_uncerts(1),3))])
            disp(['Fraction of X in bedrock fXB = ' num2str(X_MAP(end).fXB) ' ' char(177) ' ' num2str(round(X_uncerts(2),3))])  
            disp(['Fraction of calcite in bedrock fCaB = ' num2str(X_MAP(end).fCaB) ' ' char(177) ' ' num2str(round(X_uncerts(3),3))])     
    %         disp(['Fraction of quartz in bedrock fQzB = ' num2str(X_MAP.fQzB) ])
    %         disp(['Fraction of X in bedrock fXB = ' num2str(X_MAP.fXB)])  
    %         disp(['Fraction of calcite in bedrock fCaB = ' num2str(X_MAP.fCaB)])  
        case 'bedrock'
            disp(['Fraction of quartz in soil fQzS = ' num2str(X_MAP(end).fQzS) ' ' char(177) ' ' num2str(round(X_uncerts(1),3)) ])
            disp(['Fraction of X in soil fXS = ' num2str(X_MAP(end).fXS) ' ' char(177) ' ' num2str(round(X_uncerts(2),3))])
            disp(['Fraction of calcite in soil fCaS = ' num2str(X_MAP(end).fCaS) ' ' char(177) ' ' num2str(round(X_uncerts(3),3))])
    %         disp(['Fraction of quartz in soil fQzS = ' num2str(X_MAP.fQzS) ])
    %         disp(['Fraction of X in soil fXS = ' num2str(X_MAP.fXS)])
    %         disp(['Fraction of calcite in soil fCaS = ' num2str(X_MAP.fCaS)])
    end
end
% 
% export = input('Do you want to export your results? "y" or "n"? ','s');
% if strcmpi(export,'y')
%     vars = {'Name','fQzS','fCaS','fXS','fQzB','fCaB','fXB','W','D'};
%     out_table = table(txt,X.fQzS,X.fCaS,X.fXS,X.fQzB,X.fCaB,X.fXB,W,MAP/rho*10 ,'VariableNames',vars);
%     writetable(out_table,[ '.\output\JO\36Cl\' txt{1} '_' X.mode '.xlsx'])
% end
