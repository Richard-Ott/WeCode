% The script calculates the "danger zone" for a soluble target mineral. The
% "danger zone" is defined as either (1) the denudation rate below which 
% increased denudation does not lead to a decrease in nuclide concentration
% in the target mineral (D_Nmax), or (2) the denudation rate below which
% there is no unique denudation rate solution (D_unique).
% Richard Ott, 2021
%
% the current version is written for 10Be and 36Cl, could easily be
% expanded to other nuclides

clc
clear 
close all
addpath '.\subroutines'

% load data
[num,sampName,X,DEMdata] = CosmoDataRead('Test_Input_Single_Cl.xlsx');
soil_masses = 20:10:200;      % in g/cm²
W           = 10:10:100;      % soil weathering rate of soluble mineral in mm/ka
%% Calculate production rates ------------------------------------------- %

Cronus_prep = {@Cronus_prep10, @Cronus_prep36};

pars = Cronus_prep{X.n}(num,DEMdata);

%% Calculate danger zone for different soil masses and weathering rates - %
D = [10,1e3];                          % Denudation min/max in mm/ka (Dmin > Weathering rate)
thres = 0.1;                           % threshold of nuclide concentration error at which inversion has converged in % of N

D_Nmax = nan(length(W),length(soil_masses));
D_unique = nan(size(D_Nmax));
for i = 1:length(soil_masses)
    for j = 1:length(W)
        X.soil_mass = soil_masses(i);
        X.W = W(j);
        % compute D_Nmax
        [D_Nmax(j,i), ~] = solCRN_D_Nmax(pars,D,X,thres);
        try
            D_unique(j,i) = solCRN_D_unique(pars,D,X);
        catch % only possible if bedrock chemistry provided
            if X.n == 1
                D_unique(j,i) = W(j)/10*pars.sp10.rb/X.fCaB;
            elseif X.n == 2
                D_unique(j,i) = W(j)/10*pars.sp36.rb/X.fCaB;
            end
        end
    end
end

if X.n == 1
    D_Nmax = D_Nmax./pars.sp10.rb*10;
    D_unique = D_unique./pars.sp10.rb*10;
elseif X.n == 2
    D_Nmax = D_Nmax./pars.sp36.rb*10;
    D_unique = D_unique./pars.sp36.rb*10;
end % convert results to mm/ka
save('Danger_Zone_Plotdata.mat')

%% Plot the danger zone

subplot(1,2,1)
[XX,YY] = meshgrid(soil_masses,W);
cmap = load('vik.mat');
cmap = struct2cell(cmap);
cmap = cmap{1};
imagesc(soil_masses,W,D_Nmax); hold on
colormap(cmap)
contour(soil_masses,W,D_Nmax,'k','ShowText','on','LevelList',[0:10:200],'LabelSpacing',500)
set(gca,'YDir','normal');
xlabel('soil mass g/cm²');
ylabel('weathering rate mm/ka');
h = colorbar;
ylabel(h, 'Denudation rate mm/ka')

subplot(1,2,2)
imagesc(soil_masses,W,D_unique); hold on
colormap(cmap)
contour(soil_masses,W,D_unique,'k','ShowText','on','LevelList',[0:10:200],'LabelSpacing',500)
set(gca,'YDir','normal');
xlabel('soil mass g/cm²');
ylabel('weathering rate mm/ka');
h = colorbar;
ylabel(h, 'Denudation rate mm/ka')
