% This function was published within Cronus v2.1 by Marrero et al. 2016
% and addpted for WeCode by Richard Ott, 2021
%
% [erate,uncert,eratemm,uncertmm]=c14erate(sampledata,sampleuncertainties,scaling_model)
%
%  Given the data for a sample and associated one standard
%  deviation uncertainties, computes the erosion rate of the sample and uncertainty.
%
% The sampledata vector contains the following information:
%
%1. Latitude (decimal degrees, -90(S) to +90(N))
%2. Longitude (decimal degrees, 0-360 degrees east)
%3. Elevation (meters)
%4. Pressure (hPa)
%5. sample thickness (cm)
%6. bulk density (g/cm^3)
%7. Shielding factor for terrain, snow, etc. (unitless)
%8. Erosion-rate epsilon (g/(cm^2*kyr))
%9. Sample 14-C concentration (atoms of 14-C/g of target)
%10. Inheritance for 14-C (atoms 14-C/g of target)
%11. Effective attenuation length -Lambdafe (g/cm^2)
%12. Depth to top of sample (g/cm^2)
%13. Year sampled (e.g. 2010)
%
% A second input vector, sampleuncertainties, contains 1-sigma
% uncertainties for all 13 inputs.  In general, we assume that
% these 13 inputs are uncorrelated.
%
% scaling_model is one of 'DE','DU','LI','LM','SA','SF','ST' and
% informs which scaling model is being used
%
% Returns:
% erate g/(cm^2*kyr) and uncertainty.

function [erate,uncert,eratemm,uncertmm]=c14erateBasin(sampledata,sampleuncertainties,scaling_model,DEM,DB,utmzone)

%
% Make sampledata and uncertainties column vectors if they aren't already.
%
if (size(sampledata,1)==1)
  sampledata=sampledata';
end
if (size(sampleuncertainties,1)==1)
  sampleuncertainties=sampleuncertainties';
end

%
% First, check that the input data is reasonable.
%
if (length(sampledata) ~= 13)
  error('sampledata has wrong size!');
end
if (length(sampleuncertainties) ~= 13)
  error('sampleuncertainties has wrong size!');
end
%if (~exist('scaling_model','var')), scaling_model = 'all'; end
%
% Setup the physical parameters.
%
pp=physpars();
%
% Extract the sample parameters from the sampledatavector.
%
sp=samppars14(sampledata);
sp.P = stdatm(nanmean(DEM.Z(DB.Z == 1)));
%
% Get the scale factors.
%
sf=scalefacs14Basin(sp,scaling_model,DEM,DB,utmzone);
%
% We need an absolute maximum age for several purposes, including
% detecting saturated samples and setting the maximum depth for comppars.
%
maxage=50;                             % It'll be saturated after 50ka
%
% Figure out the maximum possible depth at which we'll ever need a
% production rate.  This is depthtotop + maxage * erosion (g/cm2/kyr) +
% thickness * density + a safety factor.
maxdepth=sp.depthtotop+maxage*sp.epsilon+sp.ls*sp.rb+1000;
%
% Computed parameters.
%
cp=comppars14(pp,sp,sf,maxdepth);

% Compute the nominal results. Assumes min erosion rate is 0.
%
erate=c14erateraw(pp,sp,sf,cp,scaling_model,0);
eratemm=erate/sp.rb*10;
%
%%%%%%%%
% Start off with a sum of 0.
%
uncertainty=0.0;
%
% Work through all 13 sample parameters.  For each nonzero
% uncertainty, add in the uncertainty term.
%
derivs=zeros(13,1);
for i=1:13
  if ((sampleuncertainties(i) ~= 0.0) && (~isnan(sampleuncertainties(i)))) %| (nargout > 3))
    if (sampledata(i) ~= 0.0)
      thisdelta=0.01*abs(sampledata(i));
    else
      thisdelta=0.01;
    end
    deltasampledata=sampledata;
    deltasampledata(i)=deltasampledata(i)+thisdelta;
    deltasp=samppars14(deltasampledata); deltasp.P =sp.P;
    % the next lines use the scaling factors calculated earlier in this
    % script. If you do not have uncertainty in lat,lon and elevation this
    % is fine. Otherwise use the original lines. By commenting these lines
    % and not using the catchment scaling factor calculation speeds up
    % the code.
%     deltasf=scalefacs14(deltasp,scaling_model);
    deltacp=comppars14(pp,deltasp,sf,maxdepth);
%     deltacp=comppars14(pp,deltasp,deltasf,maxdepth);
    deltaoutput=c14erateraw(pp,deltasp,sf,deltacp,scaling_model,0);
%     deltaoutput=c14erateraw(pp,deltasp,deltasf,deltacp,scaling_model,0);
    deltaerate=deltaoutput(1);
    deratei=(deltaerate-erate)/(thisdelta);
    derivs(i)=deratei;
    if (~isnan(sampleuncertainties(i)))
      uncertainty=uncertainty+(deratei^2*sampleuncertainties(i)^2);
    end
  end
end
%
% Add in terms for the uncertainty in production rates.
%
%
% Uncertainty in PsC.
%
deltapp=pp;
deltapp.PsC=pp.PsC+0.01*abs(pp.PsC);
deltaoutput=c14erateraw(deltapp,sp,sf,cp,scaling_model,0);
deltaerate=deltaoutput(1);
deratepsC=(deltaerate-erate)/(0.01*abs(pp.PsC));
uncertainty=uncertainty+(deratepsC^2*pp.sigmaPsC^2);
%
% Finally, take the square root of uncertainty to get a standard deviation.
%
uncert=sqrt(uncertainty);
%convert to other units
uncertmm=uncert/erate*eratemm;
