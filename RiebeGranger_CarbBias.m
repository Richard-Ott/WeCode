% This script calculates the bias in measurements of nuclides of minerals
% that are either more or less soluble than their host rock during soil
% weathering.
% It also calculates biases due to weathering at the soil bedrock
% interface.
% Based on Riebe and Granger eqn. 14.
% Richard Ott, 2021

% currently this cannot handle production from thermal and epithermal
% neutrons for 36Cl
clc
clear 
close all

% User choice ----------------------------------------------------------- %
nuclide = '10Be';  % nuclide of interest '10Be', '36Cl'


% PRODUCTION RATES ------------------------------------------------------ %
% attenuation lengths in g/cm²
Ln = 160; 
Lsm = 1500;
Lfm = 4320;
Pn = 4.01;
Psm = 0.1;
Pfm = 0.01;


Xs_Xr = 0:0.01:5;             % ratio of enrichment of depletion of target mineral in the soil vs bedrock Xsoil/Xrock
ph = 20:1:200;                % soil mass in g/cm²

CEF = nan(length(ph),length(Xs_Xr));  % chemical erosion factor 
p_err = CEF;                                 % percent error compared to standard denudation

% loop through all combinations of soil mass and enrichment/depletion
% ratios
for i = 1:length(ph)
    for j = 1:length(Xs_Xr)
        XX = Xs_Xr(j);
        
        % spallation for corrected and uncorrected case
        NueN = Pn*Ln.*(XX.*(1-exp(-ph(i)./Ln))+exp(-ph(i)./Ln));
        NueD = Pn*Ln;
        
        % slow muons A
        MsmN = Psm*Lsm.*(XX.*(1-exp(-ph(i)./Lsm))+exp(-ph(i)./Lsm));
        MsmD = Psm*Lsm;
        
        % fast muons
        MfmN = Pfm*Lfm.*(XX.*(1-exp(-ph(i)./Lfm))+exp(-ph(i)./Lfm));
        MfmD = Pfm*Lfm;
        
        % total N
        S_N =  NueN + MsmN + MfmN;
        S_D =  NueD + MsmD + MfmD;
        
        CEF(i,j) = S_N/S_D;
        p_err(i,j) = (CEF(i,j)-1)*100;
    end
end

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

%% Soil bedrock interface weathering

Wr = linspace(0,0.99,100);   % percentage of denudation that is weathering at soil bedrock interface

SBW = nan(length(ph),length(Wr));  % SBW bias
p_err = SBW;                       % percent error compared to standard denudation

for i = 1:length(ph)
    for j = 1:length(Wr)
        
        DrEr = 1/(1-Wr(j));
               % spallation for corrected and uncorrected case
        NueN = Pn*Ln*(DrEr-(DrEr-1)*exp(-ph(i)/Ln));
        NueD = Pn*Ln;
        
        % slow muons A
        MsmN = Psm*Lsm*(DrEr-(DrEr-1)*exp(-ph(i)/Lsm));
        MsmD = Psm*Lsm;
        
        % fast muons
        MfmN = Pfm*Lfm*(DrEr-(DrEr-1)*exp(-ph(i)/Lfm));
        MfmD = Pfm*Lfm;
        
        % total N
        S_N =  NueN + MsmN + MfmN;
        S_D =  NueD + MsmD + MfmD;
        
        SBW(i,j) = S_N/S_D;
        p_err(i,j) = (SBW(i,j)-1)*100;
        
        
    end
end



[X,Y] = meshgrid(Wr,ph);

figure()
% load colormap
cmap = load('lapaz.mat');
cmap = struct2cell(cmap);
cmap = cmap{1};
imagesc(Wr,ph,p_err,[0,500]); hold on
colormap(cmap)
levels = [0:20:200,200:100:5e3];
contour(X,Y,p_err,'LineColor','k','LevelList',levels)
contour(X,Y,p_err,[200,200],'LineColor',[.2,.2,.2],'ContourZLevel',1e6,'LineWidth',2)
set(gca,'YDir','normal');
xlabel('fraction soil-bedrock interface weathering');
ylabel('Soil mass (g/cm^2)');
h = colorbar;
ylabel(h, 'Percentage error')
ylim([min(ph),max(ph)])
xlim([min(Wr),max(Wr)])








