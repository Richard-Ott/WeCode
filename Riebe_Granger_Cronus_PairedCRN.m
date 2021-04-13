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
    num36 = num36(ind,:); txt36 = txt36(ind,:);
end

fQzB = num10(:,35); 
fCaB = num36(:,82); 

% method = 'composition_4'; % Which method do you want to use to correct your data?
% % 'composition_3' - enrichment/depletion of minerals in the soil is defined by
% % bedrock mineral composition and the concentration of one of the two
% % target minerals within the soil (3 values)
% % 'composition_4' - you have concentraions of the target minerals in both
% % bedrock and soil
% % 'weathering'    - you know an independently derived weathering rate
% % within the catchment and assume only one target mineral is soluble
% 
% % load compositional/weathering data ------------------------------------ %
% switch method
%     case 'composition_3'
% 
%     case 'composition_4'
%         fQzS = num10(:,36); 
%         fCaS = num36(:,83); 
%         fQzB = num10(:,35); 
%         fCaB = num36(:,82); 
%     case 'weathering'
% end


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
    [nominal10,uncerts10,sp10,sf10,cp10,maxage10,maxdepth10,erate_raw10] = Cronus_prep10(num,...
        scaling_model,pp,DEMdata,DEM,DB,utmzone);
    [nominal36,uncerts36,sp36,sf36,cp36,maxage36,maxdepth36,erate_raw36] = Cronus_prep36(num,...
        scaling_model,pp,DEMdata,DEM,DB,utmzone);
else
    [nominal10,uncerts10,sp10,sf10,cp10,maxage10,maxdepth10,erate_raw10] = Cronus_prep10(num,...
        scaling_model,pp,DEMdata);
    [nominal36,uncerts36,sp36,sf36,cp36,maxage36,maxdepth36,erate_raw36] = Cronus_prep36(num,...
        scaling_model,pp,DEMdata);
end

%% Compute nuclide concentrations for different enrichment/depletions and match measured values ------------- %

% ------------------------------------------------------------------- %
% Part 1: Calculate nuclid concentrations at soil-bedrock interface and
% mean production rate within soil
sp10.depthtotop = soil_mass;            % set depth to soil bedrock interface
sp36.depthtotop = soil_mass;            % set depth to soil bedrock interface
soil_depths = 1:0.1:soil_mass; 
% Figure out the maximum possible depth at which we'll ever need a
% production rate.  
maxage=500;                       % 2000ka should be saturated for 36Cl. By trial and error I found that 500 does only ~0.5% diff in the result but doubles the speed

% ADD SOME KIND OF NA-SEARCH ALGORITHM THAT HELPS FIND THE LOW MISFIT ZONE
% FAST

sp10.epsilon = erate_raw10;         % conventional erosion rate from erateraw function
sp36.epsilon = erate_raw36;         % conventional erosion rate from erateraw function
Pz10 = nan(length(soil_depths),1);  % set up empty vector for production rates in soil with depth
Pz36 = nan(length(soil_depths),1);  % set up empty vector for production rates in soil with depth

N_SBI10 = predN1026(pp,sp10,sf10,cp10,maxage,scaling_model,1);  % 10Be concentration at soil-bedrock interface 
N_SBI36 = predN36  (pp,sp36,sf36,cp36,maxage,scaling_model,1);  % 36Cl concentration at soil-bedrock interface
    
% Calculate average production rate within soil
sf10.currentsf=getcurrentsf(sf10,0,scaling_model,'be');
sf36.currentsf=getcurrentsf(sf36,0,scaling_model,'cl');
for j = 1:length(soil_depths)
    Pz10(j) = prodz1026(soil_depths(j),pp,sf10,cp10);
    Pz36(j) = prodz36  (soil_depths(j),pp,sf36,cp36);
end
P_avg10 = mean(Pz10);   % 10Be average soil production rate   
P_avg36 = mean(Pz36);   % 36Cl average soil production rate 
    
% ------------------------------------------------------------------- %
% Part 2: calculate final average soil nuclide concentrations

% final average soil mineral soil concentration
Ntot10 = N_SBI10 + P_avg10 * (soil_mass/(erate_raw10/1000)) * fQzS/fQzB;  
Ntot36 = N_SBI36 + P_avg36 * (soil_mass/(erate_raw36/1000)) * fCaS/fCaB;  

%  calculate conventional erates and see if they agree within <1%
blabliabla

CEF(j) = Ntot/num(9);    % mineral concentration/measured concentration
p_err(j) = (CEF(j)-1)*100;


%% PLOT N-ERROR AGAINST DEPLETION/ENRICHMENT ---------------------------- %
[X,Y] = meshgrid(Xs_Xr,ph);

figure()
% load colormap
cmap = load('vik.mat');
cmap = struct2cell(cmap);
cmap = cmap{1};
imagesc(Xs_Xr,ph,p_err,[-170,170]); hold on
colormap(cmap)
% contour(X,Y,p_err,'k','LevelList',[-10:5:100])
contour(X,Y,p_err,'k','ShowText','on','LevelList',[-50:10:200],'LabelSpacing',500)
set(gca,'YDir','normal');
xlabel('X_{soil}/X_{rock}');
ylabel('Soil mass (g/cm^2)');
h = colorbar;
ylabel(h, 'Percentage error')
ylim([min(ph),max(ph)])

% indicate zones of depletion and enrichment of target mineral
h = vline(1,'k-');
h.LineWidth = 2;
text(0.2,100,'Depletion','FontSize',14)
text(2.5,100,'Enrichment','FontSize',14)

