% [erate,uncert,eratemm,uncertmm]=cl36erate(sampledata,sampleuncertainties,scaling_model,DEM,DB)
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
%9. Sample 36Cl concentration (atoms of 36Cl/g of target)
%10. Inheritance for 36Cl (atoms 36Cl/g of target)
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
% DEM is a digital elevation model fully including the sample catchments as
% GRIDobj
% DB is a GRIDobj individually identifying the drainage basins of the
% measured samples as generated by the drainagebasins function.
%
% Returns:
% erate g/(cm^2*kyr) and uncertainty.

function [erate,uncert,eratemm,uncertmm]=cl36erateBasin(sampledata,sampleuncertainties,scaling_model,DEM,DB,utmzone)
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
if (length(sampledata) ~= 38)
  error('sampledata has wrong size!');
end
if (length(sampleuncertainties) ~= 38)
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
sp=samppars36(sampledata);
sp.P = stdatm(nanmean(DEM.Z(DB.Z == 1)));
%
% Get the scale factors.
%
sf=scalefacs36Basin(sp,scaling_model,DEM,DB,utmzone);
% if export
%     save(['scaling_' sampname '.mat'],'sf')
% end
%
% We need an absolute maximum age for several purposes, including
% detecting saturated samples and setting the maximum depth for comppars.
%
maxage=2000;                             % It'll be saturated after 2Ma
%
% Figure out the maximum possible depth at which we'll ever need a
% production rate.  This is depthtotop + maxage * erosion (g/cm2/kyr) +
% thickness * density + a safety factor. 
maxdepth=sp.depthtotop+maxage*sp.epsilon+sp.ls*sp.rb+1000; 
%
% Computed parameters.
%
cp=comppars36(pp,sp,sf,maxdepth);

% Compute the nominal results. Assumes min erosion rate is 0.
%
erate=cl36erateraw(pp,sp,sf,cp,scaling_model,0);

eratemm=erate/sp.rb*10;
%%%%%%%%
% Start off with a sum of 0.
%
uncertainty=0.0;
%
% Work through all 13 sample parameters.  For each nonzero
% uncertainty, add in the uncertainty term.
%
derivs=zeros(38,1);
for i=1:38
  if ((sampleuncertainties(i) ~= 0.0))&& (~isnan(sampleuncertainties(i)))
    if (sampledata(i) ~= 0.0)
      thisdelta=0.01*abs(sampledata(i));
    else
      thisdelta=0.01;
    end
    deltasampledata=sampledata;
    deltasampledata(i)=deltasampledata(i)+thisdelta;
    deltasp=samppars36(deltasampledata);
    % the next lines use the scaling factors calculated earlier in this
    % script. If you do not have uncertainty in lat,lon and elevation this
    % is fine. Otherwise use the original lines. By commenting these lines
    % and not using the catchment scaling factor calculation speeds up
    % the code.
%     deltasp.P = stdatm(nanmean(DEM.Z(DB.Z == 1)));
%     deltasf=scalefacs36Basin(deltasp,scaling_model,DEM,DB,utmzone);
%     deltacp=comppars36(pp,deltasp,deltasf,maxdepth);
    deltacp=comppars36(pp,sp,sf,maxdepth);
    deltaoutput=cl36erateraw(pp,deltasp,sf,deltacp,scaling_model,0);
%     deltaoutput=cl36erateraw(pp,deltasp,deltasf,deltacp,scaling_model,0);
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
% First, for Pf0.
%
deltapp=pp;
deltapp.Pf0=pp.Pf0+0.01*abs(pp.Pf0);
deltaoutput=cl36erateraw(deltapp,sp,sf,cp,scaling_model,0);
deltaerate=deltaoutput;
dagedpf0=(deltaerate-erate)/(0.01*abs(pp.Pf0));
uncertainty=uncertainty+(dagedpf0^2*pp.sigmapf0^2);
%
% Next, for PsCa0.
%
deltapp=pp;
deltapp.PsCa0=pp.PsCa0+0.01*abs(pp.PsCa0);
deltaoutput=cl36erateraw(deltapp,sp,sf,cp,scaling_model,0);
deltaerate=deltaoutput;
dagedpsCa0=(deltaerate-erate)/(0.01*abs(pp.PsCa0));
uncertainty=uncertainty+(dagedpsCa0^2*pp.sigmaPsCa0^2);
%
% Finally, for PsK0.
%
deltapp=pp;
deltapp.PsK0=pp.PsK0+0.01*abs(pp.PsK0);
deltaoutput=cl36erateraw(deltapp,sp,sf,cp,scaling_model,0);
deltaerate=deltaoutput;
dagedpsK0=(deltaerate-erate)/(0.01*abs(pp.PsK0));
uncertainty=uncertainty+(dagedpsK0^2*pp.sigmaPsK0^2);
%
% Finally, take the square root of uncertainty to get a standard deviation.
%
uncert=sqrt(uncertainty);
%convert to other units
uncertmm=uncert/erate*eratemm;
toc;