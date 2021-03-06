function [Ntot10, Ntot36, D, Ws, Xcur] = paired_N_forward(pp,sp10,sp36,sf10,sf36,cp10,cp36,maxage10,maxage36,scaling_model,soil_mass,x0,X)
% This functions computes the average soil nuclide concentration in a soil
% with differential weathering of 2 minerals (here assumed to be quartz and
% calcite). 
% Concentrations computed with Cronus.
% x0 = [fraction quartz, denudation rate];
% X - input composition (without the other half that needs to be computed)
% otherwise standard Cronus input for both nuclides
% Richard Ott, 2021


D = x0(2);  % denudation rate

Xcur = X;
switch X.mode
    case 'soil'
        Xcur.fQzB = x0(1); 
        Xcur.fXB  = X.fXS * (x0(1)/X.fQzS);    % the other insoluble mineral should behave like quartz
        Xcur.fCaB = 1 - x0(1) - Xcur.fXB;     % after subtracting the insoluble minerals, the rest should be the soluble mineral
    case 'bedrock'
        Xcur.fQzS = x0(1); 
        Xcur.fXS  = X.fXB * (x0(1)/X.fQzB);    % the other insoluble mineral should behave like quartz
        Xcur.fCaS = 1 - x0(1) - Xcur.fXS;     % after subtracting the insoluble minerals, the rest should be the soluble mineral
end

% set current denudation rate
sp10.epsilon = D;        
sp36.epsilon = D;        

N_SBI10 = predN1026(pp,sp10,sf10,cp10,maxage10,scaling_model,1);  % 10Be concentration at soil-bedrock interface 
N_SBI36 = predN36  (pp,sp36,sf36,cp36,maxage36,scaling_model,1);  % 36Cl concentration at soil-bedrock interface
    
% Calculate average production rate within soil
sf10.currentsf=getcurrentsf(sf10,0,scaling_model,'be');
sf36.currentsf=getcurrentsf(sf36,0,scaling_model,'cl');
% In the following scaling factors are assumed to be constant over the time
% of soil erosion. 

soil_depths = 0:0.1:soil_mass;
Pz10 = nan(1,length(soil_depths));
Pz36 = nan(1,length(soil_depths));
for j = 1:length(soil_depths)
    Pz10(j) = prodz1026(soil_depths(j),pp,sf10,cp10);
    Pz36(j) = prodz36  (soil_depths(j),pp,sf36,cp36);
end
P_avg10 = mean(Pz10);   % 10Be average soil production rate   
P_avg36 = mean(Pz36);   % 36Cl average soil production rate 
    
% ------------------------------------------------------------------- %
% Part 2: calculate final average soil nuclide concentrations

% final average soil mineral soil concentration
Ntot10 = (N_SBI10 - P_avg10/pp.lambda10Be)*exp(-pp.lambda10Be * (soil_mass/(D/1000)) * Xcur.fQzS/Xcur.fQzB)+ P_avg10/pp.lambda10Be;  
Ntot36 = (N_SBI36 - P_avg36/pp.lambda36Cl)*exp(-pp.lambda36Cl * (soil_mass/(D/1000)) * Xcur.fCaS/Xcur.fCaB)+ P_avg36/pp.lambda36Cl;  

% below equations without decay
% Ntot10 = N_SBI10 + P_avg10 * (soil_mass/(D/1000)) * Xcur.fQzS/Xcur.fQzB;
% Ntot36 = N_SBI36 + P_avg36 * (soil_mass/(D/1000)) * Xcur.fCaS/Xcur.fCaB;  

Es = Xcur.fQzB/Xcur.fQzS * D;  % soil erosion
Ws = D-Es;                     % soil weathering
end

