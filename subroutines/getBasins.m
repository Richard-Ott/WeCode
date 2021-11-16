function varargout = getBasins(DEM,x,y,coordType)
% delineate drainage basins in DEM based on x,y location
% coordType either 'map' or 'll' which stands for lat-lon
% Output drainage basin GRIDobj as derived from drainagebasins TT function.
% Richard Ott, 2021

% compute flow of DEM
FD = FLOWobj(DEM,'preprocess','carve');
A = flowacc(FD);
S = STREAMobj(FD,'minarea',1e6,'unit','map');  % you may need to adjust the minimum area for stream initiation

% snap sampling locations 2 stream
if strcmpi(coordType,'ll')
    utmzone = input('Which UTM zone is the DEM in? E.g. 35 ');
    [x_utm,y_utm] = ll2utm(y,x,utmzone);
    coords = [x_utm,y_utm]; 
else
    coords = [x,y];
end
[IX,~] = snap2stream(A > 1e3,coords);

% get drainage basins
DB = drainagebasins(FD,IX);
imageschs(DEM,DB)                          % visually check the extracted basins (can be commented)

all_good = input('Did the snapping of stream locations produce reasonable catchments "Y" or "N"? ' , 's');
if strcmpi(all_good,'N')
    error('Check your sample input locations and DEM. Something apparently went wrong with the stream snapping')
end
close all

if strcmpi(coordType,'ll')
    varargout = {DB,utmzone};
else
    varargout = DB;
end

end

