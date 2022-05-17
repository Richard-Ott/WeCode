% This function was published within Cronus v2.1 by Marrero et al. 2016
% and addpted for WeCode by Richard Ott, 2021
%
%
%  age=be10erateraw(pp,sp,sf,cp,maxerate,minrate)
%
%  Given the data for a saturated sample, computes the
%  corresponding erosion rate.
%
% This inner routine does not handle the uncertainty calculations,
% which are done by be10erate.m.  Instead, this inner routine simply
% does the basic computation of the erosion rate.
%
% Inputs:
%    pp,sp,sf,cp       as from getpars10.
%    maxrate           Maximum erosion rate (g/(cm^2*kyr))
%    minrate           Minimum erosion rate (g/(cm^2*kyr))
%                      (optional, depfaults to 0)
%
% Returns
%
%   erate              g/(cm^2*kyr)
%
%
function erate=be10erateraw(pp,sp,sf,cp,scaling_model,minrate)
%
% Set a default minimum erosion rate of 0 if none is specified.
%
if (nargin < 6)
  minrate=0.0;
end
%
% Figure out the maximum possible depth at which we'll ever need a
% production rate.
%
maxage=1000;                       %8200 ka should be saturated for 10-Be, edited to 1000
%
sp.epsilon=minrate;
maxcon=predN1026(pp,sp,sf,cp,maxage,scaling_model,1);
%
% Make sure that the concentration is reasonable.
%
if (sp.concentration10 > maxcon)
  warning('This sample is oversaturated based on min erosion rate');
  erate=NaN;
  return;
end
%
% The main loop does bisection search to find the corresponding
% erosion rate.
%
lowerrate=minrate;
%find a max erosion rate by starting at 10, checking concentration, and
%increasing by an order of magnitude until the concentration is within the
%bounds of search.
MAXUPPERRATE=300; %this number is from trial and error; put in place to
% keep code from erroring out. Sept 2016 Shasta
upperrate=MAXUPPERRATE/1000; %divides it into 3 steps by order of magnitude
sp.epsilon=upperrate;
minconc=predN1026(pp,sp,sf,cp,maxage,scaling_model,1);

upperrate = 300;
sp.epsilon=upperrate;
maxdepth=sp.depthtotop+maxage*sp.epsilon+sp.ls*sp.rb+2000;
cp=comppars1026(pp,sp,sf,maxdepth);

% while (minconc-sp.concentration10 > 0)
%     upperrate=upperrate*10;
%
%     if upperrate > MAXUPPERRATE
%        warning(['Maximum erosion rate (' num2str(MAXUPPERRATE) 'g/cm^2) reached - check sample inputs'])
%        erate=NaN;
%        return;
%     end
%
%     sp.epsilon=upperrate;
%     %recalculate the maxdepth so that cp can be found to an appropriate depth
%     maxdepth=sp.depthtotop+maxage*sp.epsilon+sp.ls*sp.rb+2000;
%     %
%     % Computed parameters.
%     %
%     cp=comppars1026(pp,sp,sf,maxdepth);
%     minconc=predN1026(pp,sp,sf,cp,maxage,scaling_model,1);
%
% end


while (upperrate-lowerrate > 1.0e-6)
  midrate=(upperrate+lowerrate)/2;
  sp.epsilon=midrate;
  midcon=predN1026(pp,sp,sf,cp,maxage,scaling_model,1);
  if (midcon > sp.concentration10)
    lowerrate=midrate;
  else
    upperrate=midrate;
  end
end
erate=(upperrate+lowerrate)/2;
