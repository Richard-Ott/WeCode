%
% sf=scalefacs1026(sample,scaling_model)
%
function sf=scalefacs1026Basin(sp1026,scaling_model,DEM,DB,utmzone)
% The DEM is assumed to be in utm coordinates whereas the sampling
% coordinates should be in lat lon.
%  Richard Ott, 2021
%
% Check if scaling_model was specified, otherwise set it to default
%
if (~exist('scaling_model','var')), scaling_model = 'all'; end

%
% Setup the scale factors.
%
sf.ST=sp1026.ST;
sf.SLth=1;
sf.SLeth=1;
sf.P=sp1026.P;
sf.elevation=sp1026.elevation;
sf.longitude=sp1026.longitude;
sf.latitude=sp1026.latitude;

load pmag_consts

IX = find(DB.Z == 1);        % get basin indices
nIX = length(IX);

% get Lat, Lon coordinates of DEM -----------------------------
% get extent of raster in UTM
x_raster = DEM.georef.SpatialRef.XWorldLimits;
y_raster = DEM.georef.SpatialRef.YWorldLimits;

% get extent of raster in Lat Lon
[raster_lat,raster_lon]=utm2ll(x_raster,y_raster,utmzone,'wgs84');

raster_lat = linspace(raster_lat(1),raster_lat(2),DEM.size(1));
raster_lon = linspace(raster_lon(1),raster_lon(2),DEM.size(2));

raster_lat = repmat(raster_lat',DEM.size(2),1);
raster_lon = repmat(raster_lon,1,DEM.size(1));


% set up empty vectors for basin values. Depending on scaling scheme very
% differnt vectors are needed.
Out = struct;

% loop through pixels in basins

wb = waitbar(0,'Calculating basin scaling factors...');
for i = 1:nIX
    tdsfsample.lat  = raster_lat(IX(i));
    tdsfsample.long = raster_lon(IX(i));
     % Perform air pressure correction for elevation (Eq. 1 from Stone, 2000)
    tdsfsample.pressure = 1013.25*exp(((-0.03417)/6.5e-3)*(log(288.15)-log(288.15-(6.5e-3*DEM.Z(IX(i))))));
    tdsfsample.elevation= DEM.Z(IX(i));
    tdsfsample.scaling=scaling_model;
    
    tdsf=get_tdsf(tdsfsample,pmag_consts);
    
    % assign output
    if (strcmpi(scaling_model,'sa') || strcmpi(scaling_model,'sf') || strcmpi(scaling_model,'all')) 

        Out(i).scaling_sp = tdsf.SF_Sf;
        Out(i).SaRc = tdsf.Rc_Sf;

        Out(i).scaling_Be = tdsf.SF_Sa10;
        Out(i).scaling_Al = tdsf.SF_Sa26;
        Out(i).scaling_He = tdsf.SF_Sa3;
        Out(i).scaling_C  = tdsf.SF_Sa14;
        Out(i).scaling_ClCa = tdsf.SF_Sa36Ca;
        Out(i).scaling_ClK  = tdsf.SF_Sa36K;
        Out(i).scaling_ClTi = tdsf.SF_Sa36Ti;
        Out(i).scaling_ClFe = tdsf.SF_Sa36Fe;
        Out(i).scaling_eth  = tdsf.SF_Saeth; 
        Out(i).scaling_th   = tdsf.SF_Sath;
        Out(i).SaRc         = tdsf.Rc_Sa;
    end
    if (strcmpi(scaling_model,'lm') || strcmpi(scaling_model,'all'))
        Out(i).SF_Lm = tdsf.SF_Lm;
        Out(i).LmRc  = tdsf.Rc_Lm;
    end
    if (strcmpi(scaling_model,'st') || strcmpi(scaling_model,'all'))
        Out(i).SF_St = tdsf.SF_St;
    end
    if (strcmpi(scaling_model,'li') || strcmpi(scaling_model,'all'))
        Out(i).SF_Li = tdsf.SF_Li;
        Out(i).LiRc  = tdsf.Rc_Li;
    end
    if (strcmpi(scaling_model,'du') || strcmpi(scaling_model,'all'))
        Out(i).SF_Du = tdsf.SF_Du;
        Out(i).DuRc  = tdsf.Rc_Du; 
    end
    if (strcmpi(scaling_model,'de') || strcmpi(scaling_model,'all'))
        Out(i).SF_De = tdsf.SF_De;
        Out(i).DeRc  = tdsf.Rc_De;
    end
    waitbar(i/nIX,wb)
end
close(wb)

sf.tdsf = tdsf;

% make the means -------------------
if (strcmpi(scaling_model,'sa') || strcmpi(scaling_model,'sf') || strcmpi(scaling_model,'all')) 

    sf.tdsf.SF_Sf = nanmean(vertcat(Out.scaling_sp),1);
    sf.tdsf.Rc_Sf = nanmean(vertcat(Out.SaRc),1);

    sf.tdsf.SF_Sa10 =  nanmean(vertcat(Out.scaling_Be),1);
    sf.tdsf.SF_Sa26 =  nanmean(vertcat(Out.scaling_Al),1);
    sf.tdsf.SF_Sa3  =  nanmean(vertcat(Out.scaling_He),1);
    sf.tdsf.SF_Sa14 =  nanmean(vertcat(Out.scaling_C),1);
    sf.tdsf.SF_Sa36Ca = nanmean(vertcat(Out.scaling_ClCa),1);
    sf.tdsf.SF_Sa36K  = nanmean(vertcat(Out.scaling_ClK),1);
    sf.tdsf.SF_Sa36Ti = nanmean(vertcat(Out.scaling_ClTi),1);
    sf.tdsf.SF_Sa36Fe = nanmean(vertcat(Out.scaling_ClFe),1);
    sf.tdsf.SF_Saeth  = nanmean(vertcat(Out.scaling_eth),1); 
    sf.tdsf.SF_Sath   = nanmean(vertcat(Out.scaling_th),1);
    sf.tdsf.Rc_Sa     = nanmean(vertcat(Out.SaRc),1);
end
if (strcmpi(scaling_model,'lm') || strcmpi(scaling_model,'all'))
    sf.tdsf.SF_Lm = nanmean(vertcat(Out.SF_Lm),1);
    sf.tdsf.Rc_Lm = nanmean(vertcat(Out.LmRc),1);
end
if (strcmpi(scaling_model,'st') || strcmpi(scaling_model,'all'))
    sf.tdsf.SF_St = nanmean(vertcat(Out.SF_St));
end
if (strcmpi(scaling_model,'li') || strcmpi(scaling_model,'all'))
    sf.tdsf.SF_Li = nanmean(vertcat(Out.SF_Li),1);
    sf.tdsf.Rc_Li = nanmean(vertcat(Out.LiRc),1);
end
if (strcmpi(scaling_model,'du') || strcmpi(scaling_model,'all'))
    sf.tdsf.SF_Du = nanmean(vertcat(Out.SF_Du),1);
    sf.tdsf.Rc_Du = nanmean(vertcat(Out.DuRc),1); 
end
if (strcmpi(scaling_model,'de') || strcmpi(scaling_model,'all'))
    sf.tdsf.SF_De = nanmean(vertcat(Out.SF_De),1);
    sf.tdsf.Rc_De = nanmean(vertcat(Out.DeRc));
end
