function Leff = neutron_att_length_DEM(DEM,utmzone)
% This script calculates the neutron attenuation length according to CRONUS
% 2.0, Marrero (2016) and Sato (2008) with a dependency on time and 
% air-pressure.
% Magnetic time evolution parameters for neutron attenuation length 
% calculation are from CRONUS 2.0 (Marrero, 2016).
% Input:        - DEM as GRIDobj
% Richard Ott, 2020

% addpath '.\CRONUS_subroutines'

%% NEUTRON ATTENUATION LENGTH ---------------------------------------------

% Perform air pressure correction for elevation (Eq. 1 from Stone, 2000)to account for diminished cosmic ray attenuation with decreased atmospheric pressue 
pressure = 1013.25*exp(((-0.03417)/6.5e-3)*(log(288.15)-log(288.15-(6.5e-3*DEM.Z)))); 

pressure_vector = pressure(:);
Z = DEM.Z(:);

Leff = nan(length(pressure_vector),1);

load pmag_consts.mat                  % get magentic field constants from CRONUS
% prepare sample file for rigiditycutoff calculation
[sample.lat,sample.long] = utm2ll(mean(DEM.georef.SpatialRef.XWorldLimits),mean(DEM.georef.SpatialRef.YWorldLimits),utmzone,'wgs84');
sample.scaling = 'sa'; 

% it takes too long to calculate the rigiditycutoff for every pixel in the
% basin. Therefore, I simply use a mean cutoff for the basin. I did a quick
% check and it seems like it doesnt really matter.
sample.elevation = nanmean(DEM.Z(:));
sample.pressure = nanmean(pressure_vector);

% rigidity cutoff
tdsf = get_tdsf(sample,pmag_consts);
Rc = mean(tdsf.Rc_Sa);

wb = waitbar(0,'Calculating effective neutron attenuation length for basin...');
for i = 1:length(pressure_vector)
    if isnan(pressure_vector(i))
        continue
    else
        Leff(i) = rawattenuationlength(pressure_vector(i),single(Rc));
    end
    waitbar(i/length(pressure_vector),wb)
end
close(wb)

Leff = nanmean(Leff);    

end

