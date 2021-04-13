% This script calculates the bias in measurements of nuclides of minerals
% that are either more or less soluble than their host rock.
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
% surface productoin rates SLHL
switch nuclide
    case '10Be'
        Pn = 4.01;          % production spallation
        % muon production following Balco 2017
        consts = struct;           % structure to store constants for muon production
        consts.f_star=0.00191;     % Model 1A, alpha=1, neg. muon capture
        consts.Natoms = 2.006e22;  % Oxygen atoms pr gram Quartz
        consts.sigma0 = 0.280e-30; % model 1A, alpha=1, fast muon interaction
        depths = 0;                % calculate muon productoin at these depths
        pressure = 1013.25;        % assume standard atmosphere at sea level
        
        Pmu = P_mu_total(depths,pressure,consts,'no');   % compute muon production rates for all depths
        
    case '36Cl'
        
end


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
        MsaN = Psm*Lsm.*(XX.*(1-exp(-ph(i)./Lsm))+exp(-ph(i)./Lsm));
        MsaD = Psm*Lsm;
        
        % slow muons B
        MsbN = Psb*Vsb.*(XX.*(1-exp(-ph(i)./Vsb))+exp(-ph(i)./Vsb));
        MsbD = Psb*Vsb;
        
        % fast muons
        MfN = Pfm*Lfm.*(XX.*(1-exp(-ph(i)./Lfm))+exp(-ph(i)./Lfm));
        MfD = Pfm*Lfm;
        
        % total N
        S_N =  NueN + MsaN + MsbN + MfN;
        S_D =  NueD + MsaD + MsbD + MfD;
        
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

