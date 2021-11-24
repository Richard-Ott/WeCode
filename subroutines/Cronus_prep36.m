function pars36 = Cronus_prep36(num,DEMdata)
% This function computes the basic structures needed for CronucsCalc
% calculations
% This function includes the max_erate_guess which should sometimes needs
% to be tuned according to your samples.

% Input: 
%           - num, the data as dervied from the Cronus input sheet
%           - DEMdata, do you use a DEM for a pixel-based productoin rate?
%           0 or 1
%           - if DEMdata == 1, then provide DEM, DB and utmzone of the DEM
%
% Note that I decreased some of the sfatey factors for the depth
% integration compared to the original Cronus calculator. For most
% denudation rates this will be ok. But if you run into errors or weird
% output, check the max_age and max_depth parameters of the depth
% integration. 
%
% Richard Ott, 2021

global scaling_model

pp=physpars();                               % get physical parameters

% First, determine the effective neutron attenuation length following
% Marrero, 2016.
if isnan(num(11)) && strcmpi(DEMdata.method,'basin')
    Leff = neutron_att_length_DEM(DEMdata.DEM,DEMdata.utmzone);
    num(11) = Leff;
elseif isnan(num(11)) && strcmpi(DEMdata.method,'location')
    Leff = neutron_att_length(num(1),num(2),num(3));
    num(11) = Leff;
end

if strcmpi(DEMdata.method,'basin')
    DEM = DEMdata.DEM;
    DB  = DEMdata.DB;
    utmzone = DEMdata.utmzone;
end

num(end+1:79) = nan;
% add default parameters ------------------------------
if isnan(num(5)), num(5) = 0.0625; end % default sample thcikness cm
if isnan(num(6)), num(6) = 2.65;   end % default sample density g/cm³
if isnan(num(7)), num(7) = 1;      end % default shielding 
if isnan(num(8)), num(8) = 0;      end % default erosion rate
if isnan(num(12)),num(12) = 0;     end % default top to sample cm
if isnan(num(79)),num(79) = 0;     end % default covariation
% -----------------------------------------------------


[nominal36,uncerts36,~]=createage36(num);    % Get basic info about the sample ages.

sp = samppars36(nominal36);
if strcmpi('basin',DEMdata.method)
    sp.P = stdatm(nanmean(DEM.Z(DB.Z == 1)));
    sf = scalefacs36Basin(sp,scaling_model,DEM,DB,utmzone);   % scaling factors
else
    if isnan(sp.P); sp.P = stdatm(sp.elevation); end % air pressure from elevation, if value not provided
    sf = scalefacs36(sp,scaling_model);   % scaling factors
end
% We need an absolute maximum age for several purposes, including
% detecting saturated samples and setting the maximum depth for comppars.
maxage36=500;            % It'll be saturated after 2Ma, decreased safety factor from original Cronus too increase speed without introducing too much bias

% Figure out the maximum possible depth at which we'll ever need a
% production rate.  This is depthtotop + maxage * erosion (g/cm2/kyr) +
% thickness * density + a safety factor. 
max_erate_guess = 400;     % maximum guess of erosion rate in area in mm/ka
maxdepth36 = maxage36*max_erate_guess+sp.ls*sp.rb+1000; 

% Computed parameters.
cp=comppars36(pp,sp,sf,maxdepth36);

% save all the parameters
pars36.nominal36   = nominal36;
pars36.uncerts36   = uncerts36;
pars36.sp36        = sp;
pars36.sf36        = sf;
pars36.cp36        = cp;
pars36.maxage36    = maxage36;
pars36.maxdepth36  = maxdepth36;
pars36.pp          = pp;

end

