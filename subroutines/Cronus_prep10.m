function pars10 = Cronus_prep10(num,DEMdata)
% This function computes the basic structures needed for CronucsCalc
% calculations
% This function includes the max_erate_guess which should sometimes needs
% to be tuned according to your samples.

% Input: 
%           - num, the data as dervied from the Cronus input sheet
%           - DEMdata, do you use a DEM for a pixel-based productoin rate?
%           0 or 1
%           - if DEMdata == 1, then provide DEM, DB and utmzone of the DEM
% Richard Ott, 2021
global scaling_model

pp=physpars();                               % get physical parameters

% First, determine the effective neutron attenuation length following
% Marrero, 2016.
if isnan(num(13)) && strcmpi(DEMdata.method,'basin')
    Leff = neutron_att_length_DEM(DEMdata.DEM,DEMdata.utmzone);
    num(13) = Leff;
elseif isnan(num(13)) && strcmpi(DEMdata.method,'location')
    Leff = neutron_att_length(num(1),num(2),num(3));
    num(13) = Leff;
end

if strcmpi(DEMdata.method,'basin')
    DEM = DEMdata.DEM;
    DB  = DEMdata.DB;
    utmzone = DEMdata.utmzone;
end

num(end+1:32) = NaN; % fills up empty values for Cronus

% add default values if empty ---------------------------------------------
if isnan(num(5)), num(5) = 0.0625; end % default sample thcikness cm
if isnan(num(6)), num(6) = 2.65;   end % default sample density g/cm³
if isnan(num(7)), num(7) = 1;      end % default shielding 
if isnan(num(8)), num(8) = 0;      end % default erosion rate
if isnan(num(14)),num(14) = 0;     end % default top to sample cm
% -------------------------------------------------------------------------


[nominal10,uncerts10]=createage1026(num);    % Get basic info about the sample ages.

sp = samppars1026(nominal10);
if strcmpi('basin',DEMdata.method)
    sp.P = stdatm(nanmean(DEM.Z(DB.Z == 1)));
    sf = scalefacs1026Basin(sp,scaling_model,DEM,DB,utmzone);   % scaling factors
else
    sp.P = stdatm(sp.elevation);            % air pressure
    sf = scalefacs1026(sp,scaling_model);   % scaling factors
end
% We need an absolute maximum age for several purposes, including
% detecting saturated samples and setting the maximum depth for comppars.
maxage=1000;                             % 8200 in original cronus, but I trhink that's just too much detail...

% Figure out the maximum possible depth at which we'll ever need a
% production rate.  This is depthtotop + maxage * erosion (g/cm2/kyr) +
% thickness * density + a safety factor. 
max_erate_guess = 400;     % maximum guess of erosion rate in area in mm/ka
maxdepth = maxage*max_erate_guess+sp.ls*sp.rb+1000; % safety factor in original cronus is 2000

% Computed parameters.
cp = comppars1026(pp,sp,sf,maxdepth);

% the denudation rate 
% erate_raw = be10erateraw(pp,sp,sf,cp,scaling_model,0);
%         eratemm=erate_raw/sp.rb*10;

pars10.nominal10   = nominal10;
pars10.uncerts10   = uncerts10;
pars10.sp10        = sp;
pars10.sf10        = sf;
pars10.cp10        = cp;
pars10.maxage      = maxage;
pars10.maxdepth    = maxdepth;
% pars10.erate_raw10 = erate_raw;
pars10.pp          = pp;

end

