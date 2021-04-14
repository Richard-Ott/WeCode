% This script calculates the denudation rate from a .
% Loosely based on Riebe and Granger eqn. 14.
% Conversly to Riebe and Granger, 2014, the nuclide concentrations are not
% calculated with exponentials. Nuclide concentraions are calculated used
% CronusCalc (Marrero, 2016). 
% The script calculates the deviations for differnent soil thicknesses
% (masses) and weathering fractions.
% Richard Ott, 2021

% the current version is written for 10Be and 36Cl, could easily be
% expanded

clc
clear 
close all
addpath 'C:\Users\r_ott\Dropbox\Richard\Crete\Cretan_fans\data'
addpath 'C:\Users\r_ott\Dropbox\Richard\PhD_ETH\matlab\CRONUS cosmo calculation\cronusearth-2.0'
addpath 'C:\Users\r_ott\Dropbox\Richard\PhD_ETH\matlab\InversionBasics\Neigborhood_Algorithm'
addpath '.\subroutines'

% User choice ----------------------------------------------------------- %
scaling_model = 'st';  % scaling model, for nomenclature see CronusCalc
[num10,txt10,~] = xlsread('10Be_data_CRONUS.xlsx',2);     % load 10Be data
[num36,txt36,~] = xlsread('36Cl_data_CRONUS.xlsx',2);     % load 36Cl data
soil_mass       = 50;                                     % average soil mass in g/cm²

DEMdata = 'location';     % Do you want to compute the erosion rate for a specific 'location', or a 'basin'
ind = input('Which of the samples would you like to run (must be the same index in both input tables?) ');
if ind ~= 0
    num10 = num10(ind,:); txt10 = txt10(ind,:);
    num36 = num36(ind+2,:); txt36 = txt36(ind+2,:);
end

X.fQzB = num10(:,35); 
X.fCaB = num36(:,82); 


%% assign data and initial basin calculations --------------------------- %

% in case your denudation rate is from an alluvial sample and you desire a
% pixel-by-pixel calculated production rate provide a DEM
if strcmpi('basin',DEMdata)
    DEM = GRIDobj();        % interactively choose the DEM that encompasses the basin
    export = 1;             % do you want to save the individual sample scaling factors as .mat file?
    % This can be useful when the computation for scaling schemes like 'sa'
    % and 'sf'  takes a long time for a big basin and you want the scaling
    % factors saved for later
    
    [DB,utmzone] = getBasins(DEM,num10(:,2),num10(:,1),'ll');  % delineate drainage basins and check their geometry
end

pp=physpars();                               % get physical parameters
nsamples=size(num10,1);                      % number of samples

%% Calculate production rates ------------------------------------------- %

if strcmpi('basin',DEMdata)
    [nominal10,uncerts10,sp10,sf10,cp10,maxage10,maxdepth10,erate_raw10] = Cronus_prep10(num10,...
        scaling_model,pp,DEMdata,DEM,DB,utmzone);
    [nominal36,uncerts36,sp36,sf36,cp36,maxage36,maxdepth36,erate_raw36] = Cronus_prep36(num36,...
        scaling_model,pp,DEMdata,DEM,DB,utmzone);
else
    [nominal10,uncerts10,sp10,sf10,cp10,maxage10,maxdepth10,erate_raw10] = Cronus_prep10(num10,...
        scaling_model,pp,DEMdata);
    [nominal36,uncerts36,sp36,sf36,cp36,maxage36,maxdepth36,erate_raw36] = Cronus_prep36(num36,...
        scaling_model,pp,DEMdata);
end

%% Compute nuclide concentrations for different enrichment/depletions and match measured values ------------- %

% ------------------------------------------------------------------- %
% set some constants
sp10.depthtotop = soil_mass;            % set depth to soil bedrock interface
sp36.depthtotop = soil_mass;            % set depth to soil bedrock interface
soil_depths = 1:0.1:soil_mass; 
% Figure out the maximum possible depth at which we'll ever need a
% production rate.  
maxage=500;                       % 2000ka should be saturated for 36Cl. By trial and error I found that 500 does only ~0.5% diff in the result but doubles the speed

% set inversion parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
itmax=200;    % Number of inversion steps
nsample=20;     % Number of samples to run every inversion (ns - in Sambridge notation, works well with ~ twice the number of dimensions)
ncells=5;      % Number of voronoi cells to re-sample, nsample-ncells will be Monte Carlo samples, (nr - in Sambridge notation, optimal range 2 - ns/2)
inv_method=0;  % Inversion method: 1 (only resampling of existing ensemble), 0 (resampling and additional random samples are generted)

% PRIORS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pnames = {'fQzS','D'};
fQzS = [0,1];                          % fraction of quartz in soil
D = [5,1e3];                           % Denudation min/max in mm/ka

D = D/10*sp10.rb;                      % convert to g/cm²/ka for Cronus, I HOPE THIS CONVERSION IS CORRECT

range_in = [fQzS',D'];                 % ranges of parameters

% Resolution, under which parameters resolution does stop the model
res_fQzS = 0.01;      
res_D = 1;                             % in mm/ka
resolution = [res_fQzS,res_D/10*sp10.rb]; 
err_max = [200;200];

nd = length(resolution);               % number of dimensions

% INIITAL MODEL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
range(:,1:nd) = range_in(:,1:nd);
scales(1:nd) = -1; % 0: No transform (All a priori model co-variances are equal to unity); 1: Use parameter range as a priori model co-variances
% -1 == Use parameter range as a priori model co-variances, see NA
% documentation
check = 0;
calcmovement=0;

% load boundaries of input parameters
load na_param_inc.mat     

% inversion parameters
monte=0;     % ???
istype=0;    %  type of initial ensemble of models used, see NA-documentation, 0 == uniform random initial sample
nsleep=1;    % ?????
noforward=0; % ?????
nclean=5000;  % ?????

% Initialize NA routines.
[restartNA,ranget,xcur]=NA_initialize(nd,nd_max,range,scales,nsample,ncells);

% Generate or read in starting models
[new_models,new_misfit] = NA_initial_sample(istype,monte,nsample,nd,range,scales);
% sort time steps
% new_models(1:steps-1,:)=sort(new_models(1:steps-1,:),'descend');  % why do the models need to be sorted???

% MAIN LOOP OF INVERSION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
misfit=[];
misfit_all=[];
models_all=[];
new_misfit=zeros(1,nsample);
obs_err=zeros(2,nsample);
models=[];
mfitmean = [];

it = 0; 
d = inf(1,nd);
% for it = 1:itmax+1
% while sum(d<resolution)<1   % run this loop until the new models achieved the resolution for both parameters
while ~any(con > 1)
    it = it+1;        
    new_misfit(:) = 0;
    obs_err(:) = 0;
    
    
    % runt the forward models for the next ns samples
    for i = 1: nsample
        X.fQzS = new_models(1,i); % fraction of quartz in soil
        X.fCaS = 1-X.fQzS;        % fraction of calcite in soil (assuming that soil only consists of calcite and quartz)
        D = new_models(2,i);    % get sample parameters
        % Run initial forward model -----------------------------------------------
        [N10m,N36m] = paired_N_forward(pp,sp10,sp36,sf10,sf36,cp10,cp36,maxage,scaling_model,soil_mass,D,X);
        
        % observational error
        obs_err(:,i) = [N10m,N36m]' - [nominal10(9),nominal36(1)]';
        new_misfit(i) = 1/2*sum((obs_err(i)./[uncerts10(9),uncerts36(1)]').^2);  % log likelihood
    end
    
    % add models and their misfit to posterior
    models=[models new_models];
    misfit=[misfit new_misfit];
    
    %	Calculate properties of current misfit distribution.
    %	(Mean,min,best model etc.)
    [mfitmean(it),mfitminc,mfitmin,mfitord,mopt]=NA_misfits(misfit);
    
    ntot=length(misfit);
    
    % tranform to scale
    models_sca = nan(size(new_models));
    for i=1:ntot
        models_sca(:,i)=transform2sca(models(:,i),nd,range,scales);
    end
    

    %	Call main NA routines %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if (monte)
        % Perform Monte Carlo search for comparison to NA.
        [new_models_sca]=NA_random(nd,range,nsample);
    else
        % generate a new sample using Neighbourhood algorithm (resample version)
        [new_models_sca,xcur,restartNA,work_NA1]=NA_sample_0517(models_sca,ntot,nsample,nd,nsleep,ncells,misfit,mfitord,ranget,check,xcur,calcmovement,nclean,inv_method);
        % I dont understand all the outputs
    end

    
    % transform the scaled models back to input units
    for i=1:nsample
        new_models(:,i)=transform2raw(nd,range,scales,new_models_sca(:,i));
    end
    
    % calculate difference between best ncells models and decide to stop
    % search
    max_d = max(models(:,mfitord(1:ncells))');
    min_d = min(models(:,mfitord(1:ncells))');
    d=max_d-min_d;
    
    con = sum(abs(obs_err) < err_max); 
    if rem(it,10) == 0
        disp(obs_err)
    end
%     if sum(d<resolution)>1  % until one of the parameters differences between cells is below 
%         % the reolsution level, keep adding models to the existing ones
%         misfit_all=[misfit_all misfit];
%         models_all=[models_all models];
%         misfit=[];
%         models=[];
%         % Generate or read in starting models
%         [new_models,new_misfit]=NA_initial_sample(istype,monte,nsample,nd,range,scales);
%         disp(num2str(it))
%     end 
end


%% OUTPUT RESULTS ------------------------------------------------------- %

figure()
% load colormap
cmap = load('vik.mat');
cmap = struct2cell(cmap);
cmap = cmap{1};


