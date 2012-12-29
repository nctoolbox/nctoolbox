function [data,grd]=nj_subsetGrid(ncRef,var,lonLatRange,dn1,dn2)
% NJ_SUBSETGRID: Subset grid based on lat-lon bounding box and time
% Usage:
%   [data, grd]=nj_subsetGrid(uri,var,[lonLatRange], [dn1], [dn2]);
% Inputs:
%   ncRef = OpenDAP Data URL, Local file, or existing 'ncgeodataset' object
%   var = variable to subset (string)
%   lonLatRange = [minLon maxLon minLat maxLat]  % matlab 'axes' function
%   dn1 = Matlab datenum, datestr or datevec
%       example: [1990 4 5 0 0 0] or '5-Apr-1990 00:00'.  If only dn1 is
%       input, the nearest time value will be returned
%   dn2 - Matlab datenum, datestr or datevec (optional).  If dn2 is also
%       input, the range of times between dn1 and dn2 will be returned
% Outputs:
%   data = subset data  based on lonlat and time range
%   grd = structure containing lon,lat and time
%
% Examples:
%   Subsetting topography/bathymetry:
%   url='http://geoport.whoi.edu/thredds/dodsC/bathy/crm_vol1.nc';
%   [data, grd]=nj_subsetGrid(url,'topo',[-71.4 -70.2 41.0 42.0]);
%
%   Subsetting model output:
%   url='http://geoport.whoi.edu/thredds/dodsC/examples/bora_feb.nc';
%   var = 'zeta';
%   lonLatRange = [13.0 16.0 41.0 42.0];  % [minlon maxlon minlat maxlat]
%   dn1 = '14-Feb-2003 12:00:00';
%   [data, grd]=nj_subsetGrid(url,var,lonLatRange,dn1); %nearest time
%   dn2 = [2003 2 16 14 0 0];
%   [data, grd]=nj_subsetGrid(url,var,lonLatRange,dn1, dn2);%time range
%
% NCTOOLBOX (http://code.google.com/p/nctoolbox)
if nargin < 2, help(mfilename), return, end
if (isa(ncRef, 'ncgeodataset')) %check for ncgeodataset Object
    nc = ncRef;
else
    nc = ncgeodataset(ncRef);
end
% get the geovariable object
geoGridVar = nc.geovariable(var);
if (~isa(geoGridVar, 'ncgeovariable'))
    disp(sprintf('MATLAB:nc_subsetGrid:Variable "%s" is not a geogrid variable.', var));
    return;
end
switch nargin
    case 2
        % read all the data (all lon/lat, all time).
        data = squeeze(geoGridVar.data);
        grd_temp = geoGridVar.grid;
    case 3
        if length(lonLatRange)==4
            structure.lat = lonLatRange([3 4]);
            structure.lon = lonLatRange([1 2]);
        elseif length(lonLatRange)==2
            structure.lat=lonLatRange(2);
            structure.lon=lonLatRange(1);
        end
        subs = geoGridVar.geosubset(structure);
        data = squeeze(subs.data);
        grd_temp = subs.grid;
    case 4
        if length(lonLatRange)==4
            structure.lat = lonLatRange([3 4]);
            structure.lon = lonLatRange([1 2]);
        elseif length(lonLatRange)==2
            structure.lat=lonLatRange(2);
            structure.lon=lonLatRange(1);
        end
        structure.time = [datenum(dn1)];
        subs = geoGridVar.geosubset(structure);
        data = squeeze(subs.data);
        grd_temp = subs.grid;
    case 5
        if length(lonLatRange)==4
            structure.lat = lonLatRange([3 4]);
            structure.lon = lonLatRange([1 2]);
        elseif length(lonLatRange)==2
            structure.lat=lonLatRange(2);
            structure.lon=lonLatRange(1);
        end
        structure.time = [datenum(dn1) datenum(dn2)];
        subs = geoGridVar.geosubset(structure);
        data = squeeze(subs.data);
        grd_temp = subs.grid;
    otherwise, error('MATLAB:nj_subsetGrid:Nargin',...
            'Incorrect number of input arguments');
end

try
    grd.lat = grd_temp.lat;
end
try
    grd.lon = grd_temp.lon;
end
try
    grd.z = grd_temp.z;
end
try
    grd.time = grd_temp.time;
end

