% This script plots figs. 3, 5, and 7 from the Ott et al., 2022
% All values computed assuming exponential production profiles with
% production rates at SLHL and Stone scaling.
% For 36Cl a pure carbonate composition is assumed and therefore only Ca
% spallation and muon production are considered.
% Richard Ott, 2021
clc; clear; close all
addpath '.\subroutines'
addpath '.\subroutines\Cronus_adaptations'

%% PRODUCTION RATES ----------------------------------------------------- %
p = 1013.25;         % hPa, standard atmosphere
Ls  = 160;           % attenuation length spallation in g/cm�
Ps10  = 4.01;        % 10Be spallation, Phillips et al., 2016
Ps36  = 52.16;       % 36Cl spallation from Ca, Marrero et al., 2016
% 10Be Muons --------
mindepth = 0; maxdepth = 7800; % g/cm2
f10_star = 0.00191;   % Model 1A, alpha=1, Balco 2017
Natoms   = 2.006e22;  % Oxygen atoms pr gram Quartz
sigma0Be = 0.280e-30; % model 1A, alpha=1, Balco 2017

p_muons=p_rate_calc2(f10_star,Natoms,sigma0Be,p,mindepth,maxdepth);

L10_fm=p_muons.L(1); % attenuation length, fast muons g/cm2
L10_nm=p_muons.L(2); % attenuation length, negative muons g/cm2
P10_fm=p_muons.P(1); % surface production rate atoms/g
P10_nm=p_muons.P(2); % surface production rate at/g

% 36Cl Muons ---------
f36_star = 0.0133582;     % Marrero et al. 2016
Natoms = 6.0169e21;       % Ca atoms per g pure CaCO3
sigma0Cl = 8.2331016e-30; % Marrero et al., 2016

p_muons=p_rate_calc2(f36_star,Natoms,sigma0Cl,p,mindepth,maxdepth);

L36_fm=p_muons.L(1); % attenuation length, fast muons g/cm2
L36_nm=p_muons.L(2); % attenuation length, negative muons g/cm2
P36_fm=p_muons.P(1); % surface production rate atoms/g
P36_nm=p_muons.P(2); % surface production rate at/g

%% SOIL WEATHERING ENRICHMENT/DEPLETION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Xs_Xr = 0:0.01:5;             % ratio of enrichment of depletion of target mineral in the soil vs bedrock Xsoil/Xrock
ph = 20:1:200;                % soil mass in g/cm�

CEF = nan(length(ph),length(Xs_Xr));      % chemical erosion factor
p_err = CEF;                              % percent error compared to standard denudation

% loop through all combinations of soil mass and enrichment/depletion
% ratios for
Ps = {Ps10, Ps36}; Pnm = {P10_nm, P36_nm}; Pfm = {P10_fm, P36_fm};
Lnm = {L10_nm, L36_nm}; Lfm = {L10_fm, L36_fm};
figure()
titles = {'10Be bias', '36Cl bias'};
% load colormap
cmap = load('vik.mat');
cmap = struct2cell(cmap);
cmap = cmap{1};
for n = 1:2
    for i = 1:length(ph)
        for j = 1:length(Xs_Xr)
            XX = Xs_Xr(j);

            % spallation for corrected and uncorrected case
            NueN = Ps{n}*Ls.*(XX.*(1-exp(-ph(i)./Ls))+exp(-ph(i)./Ls));
            NueD = Ps{n}*Ls;

            % slow muons A
            MsmN = Pnm{n}*Lnm{n}.*(XX.*(1-exp(-ph(i)./Lnm{n}))+exp(-ph(i)./Lnm{n}));
            MsmD = Pnm{n}*Lnm{n};

            % fast muons
            MfmN = Pfm{n}*Lfm{n}.*(XX.*(1-exp(-ph(i)./Lfm{n}))+exp(-ph(i)./Lfm{n}));
            MfmD = Pfm{n}*Lfm{n};

            % total N
            S_N =  NueN + MsmN + MfmN;
            S_D =  NueD + MsmD + MfmD;

            CEF(i,j) = S_N/S_D;
            p_err(i,j) = (CEF(i,j)-1)*100;
        end
    end

    [X,Y] = meshgrid(Xs_Xr,ph);

    subplot(1,2,n)
    imagesc(Xs_Xr,ph,p_err,[-170,170]); hold on
    colormap(cmap)
    % contour(X,Y,p_err,'k','LevelList',[-10:5:100])
    contour(X,Y,p_err,'k','ShowText','on','LevelList',[-50:10:250],'LabelSpacing',500)
    set(gca,'YDir','normal');
    xlabel('X_{soil}/X_{bedrock}');
    ylabel('Soil mass (g/cm^2)');
    h = colorbar;
    ylabel(h, 'Percentage error')
    ylim([min(ph),max(ph)])
    title(titles{n})

    % indicate zones of depletion and enrichment of target mineral
    h = vline(1,'k-');
    h.LineWidth = 2;
    text(0.2,100,'Depletion','FontSize',14)
    text(2.5,100,'Enrichment','FontSize',14)
end
%% Soil bedrock interface weathering %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Wr = linspace(0,0.99,100);   % percentage of denudation that is weathering at soil bedrock interface

SBW = nan(length(ph),length(Wr));  % SBW bias
p_err = SBW;                       % percent error compared to standard denudation

figure()
for n = 1:2
    for i = 1:length(ph)
        for j = 1:length(Wr)

            DrEr = 1/(1-Wr(j));
                   % spallation for corrected and uncorrected case
            NueN = Ps{n}*Ls*(DrEr-(DrEr-1)*exp(-ph(i)/Ls));
            NueD = Ps{n}*Ls;

            % slow muons A
            MsmN = Pnm{n}*Lnm{n}*(DrEr-(DrEr-1)*exp(-ph(i)/Lnm{n}));
            MsmD = Pnm{n}*Lnm{n};

            % fast muons
            MfmN = Pfm{n}*Lfm{n}*(DrEr-(DrEr-1)*exp(-ph(i)/Lfm{n}));
            MfmD = Pfm{n}*Lfm{n};

            % total N
            S_N =  NueN + MsmN + MfmN;
            S_D =  NueD + MsmD + MfmD;

            SBW(i,j) = S_N/S_D;
            p_err(i,j) = (SBW(i,j)-1)*100;


        end
    end

    subplot(1,2,n)
    [X,Y] = meshgrid(Wr,ph);
    % load colormap
    cmap = load('lapaz.mat');
    cmap = struct2cell(cmap);
    cmap = cmap{1};
    imagesc(Wr,ph,p_err,[0,500]); hold on
    colormap(cmap)
    levels = [0:20:200,200:100:5e3];
    contour(X,Y,p_err,'LineColor','k','LevelList',levels,'LineWidth',1)
    contour(X,Y,p_err,[200,200],'LineColor',[.2,.2,.2],'ContourZLevel',1e6,'LineWidth',2.2)
    set(gca,'YDir','normal');
    xlabel('fraction soil-bedrock interface weathering');
    ylabel('Soil mass (g/cm^2)');
    h = colorbar;
    ylabel(h, 'Percentage error')
    ylim([min(ph),max(ph)])
    xlim([min(Wr),max(Wr)])
    title(titles{n})
end


%% Danger zone plot for pure carbonate at SLHL %%%%%%%%%%%%%%%%%%%%%%%%%%%%

Xb = 0.7;
figure()
% plot the N-curve for 1 sample -------------------------------------------
subplot(1,3,1)
Wsamp = 13.5; %g /cm�/ka 13.5
ph_samp = 200; % 200
D = 0:1:100;
D = D./1e3; Wsamp = Wsamp./1e3; % convert to g/cm/a
Ns  = Ps36  .*Ls    ./D.*exp(-ph_samp/Ls)     + (Ps36   .*(1-exp(-ph_samp./Ls))     .* ph_samp .*(1-Wsamp./(D.*Xb)))./ (D-Wsamp);
Nnm = P36_nm.*L36_nm./D.*exp(-ph_samp/L36_nm) + (P36_nm .*(1-exp(-ph_samp./L36_nm)) .* ph_samp .*(1-Wsamp./(D.*Xb)))./ (D-Wsamp);
Nfm = P36_fm.*L36_fm./D.*exp(-ph_samp/L36_fm) + (P36_fm .*(1-exp(-ph_samp./L36_fm)) .* ph_samp .*(1-Wsamp./(D.*Xb)))./ (D-Wsamp);
Ntot = Ns + Nnm + Nfm;
Dmin = Wsamp/Xb;
inds = D < Dmin;
Ntot(inds) = nan;
plot(D/2.65*1e4,Ntot,'k-','LineWidth',1.5);
xlim([0, inf])
xlabel('Denudation rate mm/ka')
ylabel('nuclide concentration at/g')

% plot D_Nmax ----------------------------------------------------------- %
W = 2.65:0.2:26.5;   % g/cm2/ka
W = W./1e3;          % g/cm2/a
D_Nmax = nan(length(ph),length(W));
D_unique = size(D_Nmax);
D_unique_dmin = false(size(D_Nmax));
options = optimset('MaxIter',5e4,'TolFun',100,'TolX',0.05);
for i = 1:length(ph)
    for j = 1:length(W)

        % D_Nmax %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % this is derivative of nuclide concentrations for a soluble
        % mineral
        Ndiv = @(x) -Ps36.*Ls     .*x.^(-2) *exp(-ph(i)/Ls)    - (Ps36   .*(1-exp(-ph(i)./Ls))     .* ph(i) .*(Xb.*x.^2-2.*W(j).*x+ W(j)^2))./ (Xb.*x.^2.*(x-W(j)).^2)...
        -P36_nm*L36_nm.*x.^(-2)*exp(-ph(i)/L36_nm) - (P36_nm .*(1-exp(-ph(i)./L36_nm)) .* ph(i) .*(Xb.*x.^2-2.*W(j).*x+ W(j)^2))./ (Xb.*x.^2.*(x-W(j)).^2)...
        -P36_fm*L36_fm.*x.^(-2)*exp(-ph(i)/L36_fm) - (P36_fm .*(1-exp(-ph(i)./L36_fm)) .* ph(i) .*(Xb.*x.^2-2.*W(j).*x+ W(j)^2))./ (Xb.*x.^2.*(x-W(j)).^2);

        D_Nmax(i,j) = fzero(Ndiv,1);
        if D_Nmax(i,j)< W(j)
            disp('too low value, finding higher solution')
            D_Nmax(i,j) = fzero(Ndiv,2);
        end

        % if fzero fails, try a lower starting position...
        if isnan(D_Nmax(i,j))
            D_Nmax(i,j) = fzero(Ndiv,0.1);
        end
        if D_Nmax(i,j)< W(j)
            disp('too low value, finding higher solution')
            D_Nmax(i,j) = fzero(Ndiv,0.2);
        end

        x0 = 0.1;
        while isnan(D_Nmax(i,j))
            x0 = x0/2;
            try
                D_Nmax(i,j) = fzero(Ndiv,x0);
                while D_Nmax(i,j) < W(j)
                    x0 = x0 + 0.002;
                    D_Nmax(i,j) = fzero(Ndiv,x0);
                end
            catch
                x0 = x0+1;
                D_Nmax(i,j) = fzero(Ndiv,x0);
                while D_Nmax(i,j) < W(j)
                    x0 = x0 + 0.002;
                    D_Nmax(i,j) = fzero(Ndiv,x0);
                end
           end
        end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Dmin = W(j)/Xb;  % minimum dendudation rate
        N = @(x) Ps36  .*Ls    ./x.*exp(-ph(i)/Ls)     + (Ps36   .*(1-exp(-ph(i)./Ls))     .* ph(i) .*(1-W(j)./(x.*Xb)))./ (x-W(j))+...
                 P36_nm.*L36_nm./x.*exp(-ph(i)/L36_nm) + (P36_nm .*(1-exp(-ph(i)./L36_nm)) .* ph(i) .*(1-W(j)./(x.*Xb)))./ (x-W(j))+...
                 P36_fm.*L36_fm./x.*exp(-ph(i)/L36_fm) + (P36_fm .*(1-exp(-ph(i)./L36_fm)) .* ph(i) .*(1-W(j)./(x.*Xb)))./ (x-W(j));

        N_Dmin = N(Dmin); % nuclide concentration Dmin
        N = @(x) Ps36  .*Ls    ./x.*exp(-ph(i)/Ls)     + (Ps36   .*(1-exp(-ph(i)./Ls))     .* ph(i) .*(1-W(j)./(x.*Xb)))./ (x-W(j))+...
                 P36_nm.*L36_nm./x.*exp(-ph(i)/L36_nm) + (P36_nm .*(1-exp(-ph(i)./L36_nm)) .* ph(i) .*(1-W(j)./(x.*Xb)))./ (x-W(j))+...
                 P36_fm.*L36_fm./x.*exp(-ph(i)/L36_fm) + (P36_fm .*(1-exp(-ph(i)./L36_fm)) .* ph(i) .*(1-W(j)./(x.*Xb)))./ (x-W(j))...
                 - N_Dmin;


        try
            D_unique(i,j) = fzero(N,[Dmin+0.00001,D(end)]);

        catch
            D_unique(i,j) = Dmin;
            D_unique_dmin(i,j) = true;
        end
   end
end
% convert back to mm/ka
W = W./2.65*1e4;
D_Nmax = D_Nmax./2.65*1e4;
D_unique = D_unique./2.65*1e4;

% plot D_Nmax ---------------------
subplot(1,3,2)
[X,Y] = meshgrid(W,ph);
% load colormap
cmap = load('lapaz.mat');
cmap = struct2cell(cmap);
cmap = cmap{1};
imagesc(W,ph,D_Nmax); hold on
colormap(cmap)
levels = [0:20:200];
contour(X,Y,D_Nmax,'LineColor','k','LevelList',levels,'ShowText','on')
set(gca,'YDir','normal');
xlabel('weathering rate mm/ka');
ylabel('Soil mass (g/cm^2)');
h = colorbar;
ylabel(h, 'D_Nmax')
ylim([min(ph),max(ph)])
xlim([min(W),max(W)])

% plot D_unique ---------------------
subplot(1,3,3)
[X,Y] = meshgrid(W,ph);
% load colormap
cmap = load('lapaz.mat');
cmap = struct2cell(cmap);
cmap = cmap{1};
imagesc(W,ph,D_unique); hold on
colormap(cmap)
levels = [0:20:400];
contour(X,Y,D_unique,'LineColor','k','LevelList',levels,'ShowText','on')
set(gca,'YDir','normal');
xlabel('weathering rate mm/ka');
ylabel('Soil mass (g/cm^2)');
h = colorbar;
ylabel(h, 'D_Nmax')
ylim([min(ph),max(ph)])
xlim([min(W),max(W)])
