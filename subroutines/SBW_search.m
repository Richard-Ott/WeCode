function D = SBW_search(pars,X)
% Calculates the corrected denudation rate in case of soil-bedrock
% interface weathering
% Richard Ott, 2021

global scaling_model

% Define options for optimization  
options = optimset('MaxIter',5e4,'TolFun',5e2,'TolX',0.05);            
% These options may need to be tuned specifically to your problem
% TolFun, maximum value that the function is allowed to be off, at/g
% TolX, the tolerance value in x-direction (erosion rates in g/cm2/a), which the
% function will use for subsequent iterations as tolerance level for
% stopping the algorithm. 


switch X.n
    case 1
        fun = @(x) abs(N10_SBW_forward(pars,scaling_model,X.soil_mass,x,X.W) - pars.nominal10(9));
        x0 = be10erateraw(pars.pp,pars.sp10,pars.sf10,pars.cp10,scaling_model,0);
        [D,~,~] = fminsearch(fun,x0,options); 
        D = D/pars.sp10.rb *10;
    case 2
        fun = @(x) abs(N36_SBW_forward(pars,scaling_model,X.soil_mass,x,X.W) - pars.nominal36(1));
        x0 = cl36erateraw(pars.pp,pars.sp36,pars.sf36,pars.cp36,scaling_model,0);
        [D,~,~] = fminsearch(fun,x0,options); 
        D = D/pars.sp36.rb *10;
end



end

