function [Nm,Xcur] = Comp_and_N_forwardX(pp,sp,sf,cp,maxage,soil_mass,D,X)
% Calculates the forward models for a nuclide given a denudation rate D, a
% weathering rate X.W, a input composition X, and the standard Cronus CRN
% sample parameters. Output is the modelled nuclide concentration and final
% composition.
% Richard Ott, 2021

global scaling_model
W = X.W/10*sp.rb;                        % convert to g/cm²/ka 

% Compute composition -----------------------------------------------
Xcur  = X;
fE = 1 - W/D;                          % fraction of erosion (to total denudation)
switch X.mode
    case 'soil'
        Xcur.fQzB = Xcur.fQzS * fE;           % quartz fraction 
        Xcur.fXB  = Xcur.fXS  * fE;           % clay fraction
        Xcur.fCaB = 1 - Xcur.fQzB - Xcur.fXB;
    case 'bedrock'
        Xcur.fQzS = Xcur.fQzB * (1/fE);
        Xcur.fXS  = Xcur.fXB  * (1/fE); 
        Xcur.fCaS = 1 - Xcur.fQzS - Xcur.fXS;
end

forward_model = {@N10_forward,@N36_forward};        % to avoid opening more switch statements, I compile the functions for each case here

[Nm,~] = forward_model{X.n}(pp,sp,sf,cp,maxage,scaling_model,soil_mass,D,Xcur);

end

