function Ntot36 = N36_SBW_forward(pars,scaling_model,soil_mass,D,W)
% This functions computes the average soil nuclide concentration in a soil
% with weathering at the soil-bedrock interface.
% Concentrations computed with Cronus.
% Richard Ott, 2021

v2struct(pars)
W = W/10*sp36.rb; % convert from mm/ka to g/cm2/ka

% set current denudation rate
sp36.epsilon = D;  
sp36.depthtotop = soil_mass;           % set depth to soil bedrock interface

N_SBI36 = predN36(pp,sp36,sf36,cp36,maxage36,scaling_model,1);  % 36Cl concentration at soil-bedrock interface
    
% Calculate average production rate within soil
sf36.currentsf=getcurrentsf(sf36,0,scaling_model,'cl');
% In the following scaling factors are assumed to be constant over the time
% of soil erosion. 

soil_depths = 0:0.1:soil_mass;
Pz36 = nan(1,length(soil_depths));
for j = 1:length(soil_depths)
    Pz36(j) = prodz36  (soil_depths(j),pp,sf36,cp36);
end
P_avg36 = mean(Pz36);   % 36Cl average soil production rate 
    
% ------------------------------------------------------------------- %
% Part 2: calculate final average soil nuclide concentrations

% final average soil mineral soil concentration
Ntot36 = N_SBI36 + P_avg36 * (soil_mass/((D-W)/1000));  
end

