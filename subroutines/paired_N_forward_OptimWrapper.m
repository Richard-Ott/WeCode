function obs_err = paired_N_forward_OptimWrapper(pp,sp10,sp36,sf10,sf36,cp10,cp36,maxage10,maxage36,scaling_model,soil_mass,x0,X,Nobs)
% wrapper function for pairedCRN optimization that minimizes the misfit
% between the forward model function (paired_N_forward) and the measured
% nuclide concentrations of Nobs = [N10;N36]
% Richard Ott, 2021


[Ntot10, Ntot36, ~, ~, ~] = paired_N_forward(pp,sp10,sp36,sf10,sf36,cp10,cp36,maxage10,maxage36,scaling_model,soil_mass,x0,X);


obs_err = sum(abs(([Ntot10;Ntot36] - Nobs)./Nobs)); % error of both combined

end

