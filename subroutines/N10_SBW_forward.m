function Ntot10 = N10_SBW_forward(pars,scaling_model,soil_mass,D,W)
% This functions computes the average soil nuclide concentration in a soil
% with weathering at the soil-bedrock interface.
% Concentrations computed with Cronus.
% Richard Ott, 2021

v2struct(pars)
W = W/10*sp10.rb;       % convert from mm/ka to g/cm2/ka

% set current denudation rate
sp10.epsilon = D;  
sp10.depthtotop = soil_mass;           % set depth to soil bedrock interface

N_SBI10 = predN1026(pp,sp10,sf10,cp10,maxage10,scaling_model,1);  % 10Be concentration at soil-bedrock interface 
    
% Calculate average production rate within soil
sf10.currentsf=getcurrentsf(sf10,0,scaling_model,'be'); % this extracts modern scaling factors.
% In the following scaling factors are assumed to be constant over the time
% of soil erosion. 

soil_depths = 0:0.1:soil_mass;
Pz10 = nan(1,length(soil_depths));
for j = 1:length(soil_depths)
    Pz10(j) = prodz1026(soil_depths(j),pp,sf10,cp10);
end
P_avg10 = mean(Pz10);   % 10Be average soil production rate   
    
% ------------------------------------------------------------------- %
% Part 2: calculate final average soil nuclide concentrations

% final average soil mineral soil concentration
Ntot10 = N_SBI10 + P_avg10 * (soil_mass/((D-W)/1000));  
end

