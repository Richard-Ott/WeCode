% [erate,eratemm]=be10erate(sampledata,sampleuncertainties,scaling_model)
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
%9. Sample 10-Be concentration (atoms of 10-Be/g of target)
%10. Inheritance for 10-Be (atoms 10-Be/g of target)
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

function [erate,uncerts,eratemm,uncertsmm]=be10erateBasin(sampledata,sampleuncertainties,scaling_model,DEM,DB,utmzone)
tic;
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
if (length(sampledata) ~= 15)
  error('sampledata has wrong size!');
end
if (length(sampleuncertainties) ~= 15)
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
sp=samppars1026(sampledata);
sp.P = stdatm(nanmean(DEM.Z(DB.Z == 1)));

%
% Get the scale factors.
%
sf=scalefacs1026Basin(sp,scaling_model,DEM,DB,utmzone);
%
% We need an absolute maximum age for several purposes, including
% detecting saturated samples and setting the maximum depth for comppars.
%
maxage=1000;                             % It'll be saturated after 8200ka, edited 1000 for speed
%
% Figure out the maximum possible depth at which we'll ever need a
% production rate.  This is depthtotop + maxage * erosion (g/cm2/kyr) +
% thickness * density + a safety factor. 
maxdepth=sp.depthtotop+maxage*sp.epsilon+sp.ls*sp.rb+2000;
%
% Computed parameters.
%
cp=comppars1026(pp,sp,sf,maxdepth);

% Compute the nominal results. Assumes min erosion rate is 0.
%
erate=be10erateraw(pp,sp,sf,cp,scaling_model,0);
if isnan(erate)
    uncerts=nan;
    eratemm=nan;
    uncertsmm=nan;
    return
end

eratemm=erate/sp.rb*10;
%%%%%%%% Uncertainty calcs
% Start off with a sum of 0.
%
uncertainty=0.0;
%
% Work through all 15 sample parameters.  For each nonzero
% uncertainty, add in the uncertainty term. Should skip Al concentration
% unc.
%
derivs=zeros(15,1);
for i=1:15
  if ((sampleuncertainties(i) ~= 0.0)) && (~isnan(sampleuncertainties(i)))
    if (sampledata(i) ~= 0.0)
      thisdelta=0.01*abs(sampledata(i));
    else
      thisdelta=0.01;
    end
    deltasampledata=sampledata;
    deltasampledata(i)=deltasampledata(i)+thisdelta;
    deltasp=samppars1026(deltasampledata);
    % the next lines use the scaling factors calculated earlier in this
    % script. If you do not have uncertainty in lat,lon and elevation this
    % is fine. Otherwise use the original lines. By commenting these lines
    % and not using the catchment scaling factor calculation speeds up
    % the code.
%     deltasf=scalefacs1026Basin(deltasp,scaling_model,DEM,DB,utmzone);
    deltacp=comppars1026(pp,deltasp,sf,maxdepth);
%     deltacp=comppars1026(pp,deltasp,deltasf,maxdepth);
    deltaoutput=be10erateraw(pp,deltasp,sf,deltacp,scaling_model,0);
%     deltaoutput=be10erateraw(pp,deltasp,deltasf,deltacp,scaling_model,0);
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
% Uncertainty in PsBe.
%
deltapp=pp;
deltapp.PsBe=pp.PsBe+0.01*abs(pp.PsBe);
deltaoutput=be10erateraw(deltapp,sp,sf,cp,scaling_model,0);
deltaerate=deltaoutput(1);
deratepsBe=(deltaerate-erate)/(0.01*abs(pp.PsBe));
uncertainty=uncertainty+(deratepsBe^2*pp.sigmaPsBe^2);
%
% Finally, take the square root of uncertainty to get a standard deviation.
%
uncerts=sqrt(uncertainty);
%convert to other units
uncertsmm=uncerts/erate*eratemm;
toc;