function [Ntot36, D_CaS] = N36_forward(pp,sp36,sf36,cp36,maxage,scaling_model,soil_mass,D,X)
% This functions computes the average soil 36Cl concentration in a soil
% with differential weathering of 2 minerals (here assumed to be quartz and
% calcite). 
% Concentrations computed with Cronus.
% Richard Ott, 2021


% set current denudation rate
sp36.epsilon = D;        

N_SBI36 = predN36(pp,sp36,sf36,cp36,maxage,scaling_model,1);  % 36Cl concentration at soil-bedrock interface
    
% Calculate average production rate within soil
sf36.currentsf=getcurrentsf(sf36,0,scaling_model,'cl');

soil_depths = 0:0.1:soil_mass;
Pz36 = nan(1,length(soil_depths));
for j = 1:length(soil_depths)
    Pz36(j) = prodz36  (soil_depths(j),pp,sf36,cp36);
end
P_avg36 = mean(Pz36);   % 36Cl average soil production rate 
    
% ------------------------------------------------------------------- %
% Part 2: calculate final average soil nuclide concentrations

% final average soil mineral soil concentration
Ntot36 = N_SBI36 + P_avg36 * (soil_mass/(D/1000)) * X.fCaS/X.fCaB;  

D_CaS = D * X.fCaB/X.fCaS;

end

