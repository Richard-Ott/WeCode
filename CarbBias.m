% This script calculates the bias in measurements for minerals
% that are either more or less soluble than their host rock.
% The script calculates for different weathering fractions and soil
% thicknesses how the bias would affect this specific sample.
%
% Conversly to Riebe and Granger, 2014, the nuclide concentrations are not
% calculated with exponentials. Nuclide concentraions are calculated used
% CronusCalc (Marrero, 2016). 
% The script calculates the deviations for differnent soil thicknesses
% (masses) and weathering fractions.
% Richard Ott, 2021

% the current version is written for 10Be and 36Cl, could easily be
% expanded

clc; clear; close all
addpath '.\subroutines'
addpath '.\subroutines\Cronus_adaptations'

% User choice ----------------------------------------------------------- %
nuclide = '36Cl';      % nuclide of interest e.g. '10Be', '36Cl'
scaling_model = 'st';  % scaling model, for nomenclature see CronusCalc
[num,txt,~] = xlsread('Test_Input_Single_Cl.xlsx',1);
DEMdata.method = 'location';     % Do you want to compute the erosion rate for a specific 'location', or a 'basin'

%% assign data and initial basin calculations --------------------------- %

% in case your denudation rate is from an alluvial sample and you desire a
% pixel-by-pixel calculated production rate provide a DEM
if strcmpi('basin',DEMdata.method)
    DEMdata.DEM = GRIDobj();        % interactively choose the DEM that encompasses the basin
    % Scaling schemes like 'sa' and 'sf'  can take a long time for a big
    % basin and you want the scaling. You may want to save the scaling data
    % for later re-runs.

    [DEMdata.DB,DEMdata.utmzone] = getBasins(DEMdata.DEM,num(:,2),num(:,1),'ll');  % delineate drainage basins and check their geometry
end

pp=physpars();                             % get physical parameters
nsamples=size(num,1);                      % number of samples


%% Calculate production rates ------------------------------------------- %

switch nuclide; case '10Be'; n = 1; case '36Cl'; n = 2; end
Cronus_prep = {@Cronus_prep10, @Cronus_prep36};
pars = Cronus_prep{n}(num,DEMdata);
v2struct(pars)

% rename some variables.
if exist('sp10')
    nominal = nominal10; uncerts = uncerts10; sp = sp10; sf = sf10; cp = cp10;
    clear nominal10 uncerts10 sp10 sf10 cp10
elseif exist('sp36')
    nominal = nominal36; uncerts = uncerts36; sp = sp36; sf = sf36; cp = cp36;
    clear nominal36 uncerts36 sp36 sf36 cp36
end

erate_funcs = {@be10erateraw, @cl36erateraw};
erate_raw=erate_funcs{n}(pp,sp,sf,cp,scaling_model,0);

%% Compute for different soil masses and enrichment factors ------------- %
Xs_Xr = 0:0.01:5;             % ratio of enrichment/depletion of target mineral in the soil vs bedrock Xsoil/Xrock
ph = 20:1:200;                % soil mass in g/cm�

CEF = nan(length(ph),length(Xs_Xr));         % chemical erosion factor
p_err = CEF;                                 % percent error compared to standard denudation

% loop through all combinations of soil mass and enrichment/depletion
% ratios
for i = 1:length(ph)

    % ------------------------------------------------------------------- %
    % Part 1: Calculate nuclid concentrations at soil-bedrock interface and
    % mean production rate within soil
    sp.depthtotop = ph(i);            % set depth to soil bedrock interface
    sp.epsilon = erate_raw;           % conventional erosion rate from erateraw function
    soil_depths = 0:0.1:ph(i);
    Pz = nan(length(soil_depths),1);  % set up empty vector for production rates in soil with depth

    switch nuclide
            case '10Be'
                % Figure out the maximum possible depth at which we'll ever need a
                % production rate.
                maxage=500;                       % 2000ka should be saturated for 36Cl. By trial and error I found that 500 does only ~0.5% diff in the result but doubles the speed

                N_SBI = predN1026(pp,sp,sf,cp,maxage,scaling_model,1);  %concentration at soil-bedrock interface

                % Calculate average production rate within soil
                sf.currentsf=getcurrentsf(sf,0,scaling_model,'cl');
                for j = 1:length(soil_depths)
                    Pz(j) = prodz1026(soil_depths(j),pp,sf,cp);
                end
                P_avg = mean(Pz);   % average soil production rate

            case '36Cl'
                % Figure out the maximum possible depth at which we'll ever need a
                % production rate.
                maxage=500;                       % 2000ka should be saturated for 36Cl. By trial and error I found that 500 does only ~0.5% diff in the result but doubles the speed

                N_SBI = predN36(pp,sp,sf,cp,maxage,scaling_model,1);  %concentration at soil-bedrock interface

                % Calculate average production rate within soil
                sf.currentsf=getcurrentsf(sf,0,scaling_model,'cl');
                for j = 1:length(soil_depths)
                    Pz(j) = prodz36(soil_depths(j),pp,sf,cp);
                end
                P_avg = mean(Pz);   % average soil production rate
    end

    % ------------------------------------------------------------------- %
    % Part 2: loop through mineral enrichment/depletion fractions and
    % calculate final average soil nuclide concentrations
    for j = 1:length(Xs_Xr)
        XX = Xs_Xr(j);

        % final average soil mineral soil concentration
        Ntot = N_SBI + P_avg * (ph(i)/(erate_raw/1000)) * Xs_Xr(j);

        CEF(i,j) = Ntot/num(9);    % mineral concentration/measured concentration
        p_err(i,j) = (CEF(i,j)-1)*100;
    end
end


%% PLOT N-ERROR AGAINST DEPLETION/ENRICHMENT ---------------------------- %
[X,Y] = meshgrid(Xs_Xr,ph);

figure()
% load colormap
cmap = load('vik.mat');
cmap = struct2cell(cmap);
cmap = cmap{1};
imagesc(Xs_Xr,ph,p_err,[-170,170]); hold on
colormap(cmap)
contour(X,Y,p_err,'k','ShowText','on','LevelList',[-50:10:200],'LabelSpacing',500)
set(gca,'YDir','normal');
xlabel('X_{regolith}/X_{rock}');
ylabel('regolith mass (g/cm^2)');
h = colorbar;
ylabel(h, 'Percentage error')
ylim([min(ph),max(ph)])

% indicate zones of depletion and enrichment of target mineral
h = vline(1,'k-');
h.LineWidth = 2;
text(0.2,100,'Depletion','FontSize',14)
text(2.5,100,'Enrichment','FontSize',14)
