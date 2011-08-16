function [data,grd]=nj_subsetGrid(ncRef,var,lonLatRange,dn1,dn2)
% NJ_SUBSETGRID: Subset grid based on lat-lon bounding box and time
%
% Usage:
%   [data, grd]=nj_subsetGrid(uri,var,[lonLatRange], dn1, dn2);
% where,
%   ncRef - Reference to netcdf file. It can be either of two
%           a. local file name or a URL  or
%           b. An 'mDataset' matlab object, which is the reference to already
%              open netcdf file.
%              [ncRef=mDataset(uri)]
%
%  var - variable to subset
%  lonLatRange - [minLon maxLon minLat maxLat]  % matlab 'axes' function
%  dn1 - Matlab datenum, datestr or datevec  ex: [1990 4 5 0 0 0] or '5-Apr-1990 00:00'
%  dn2 - Matlab datenum, datestr or datevec (optional)
%
% Returns,
%   data - subset data  based on lonlat and time range
%   grd - structure containing lon,lat and time
%
%  e.g,
%   ncRef='http://www.gri.msstate.edu/rsearch_data/nopp/bora_feb.nc';
%   var = 'temp';
%   lonLatRange = [13.0 16.0 41.0 42.0];   [minlon maxlon minlat maxlat]
%   dn1 = '14-Feb-2003 12:00:00';
%   dn2 = [2003 2 16 14 0 0];
%   [data, grd]=nj_subsetGrid(ncRef,var,lonLatRange,dn1, dn2)     or
%   [data, grd]=nj_subsetGrid(ncRef,var,lonLatRange,dn1)          or
%   [data, grd]=nj_subsetGrid(ncRef,var,lonLatRange)              or
%   [data, grd]=nj_subsetGrid(ncRef,var)                          or
%
%
%
% Sachin Kumar Bhate (skbhate@ngi.msstate.edu)  (C) 2008
% Mississippi State University%

% import the NetCDF-Java methods

if nargin < 2, help(mfilename), return, end

%initialize
%data (volume or subset)
data=[];
%structure containing lon,lat
grd.lat=[];
grd.lon=[];
grd.time=[];
grd.z=[];
jdmat = [];
isNcRef=0;


if (isa(ncRef, 'ncgeodataset')) %check for mDataset Object
  nc = ncRef;
  isNcRef=1;
else
  % open CF-compliant NetCDF File as a Common Data Model (CDM) "Grid Dataset"
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
    myData = squeeze(geoGridVar.data);
    grd = geoGridVar.grid;
  case 3
    if length(lonLatRange)==4
      structure.lat = lonLatRange([3 4]);
      structure.lon = lonLatRange([1 2]);
    elseif length(lonLatRange)==2
      structure.lat=lonLatRange(2);
      structure.lon=lonLatRange(1);
    end
    subs = geoGridVar.geosubset(structure);
    myData = squeeze(subs.data);
    grd=subs.grid;
  case 4
    if length(lonLatRange)==4
      structure.lat = lonLatRange([3 4]);
      structure.lon = lonLatRange([1 2]);
    elseif length(lonLatRange)==2
      structure.lat=lonLatRange(2);
      structure.lon=lonLatRange(1);
    end
    structure.time = [datenum(dn1) datenum(dn1)];
    subs = geoGridVar.geosubset(structure);
    myData = squeeze(subs.data);
    grd=subs.grid;
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
    myData = squeeze(subs.data);
    grd=subs.grid;
    
  otherwise, error('MATLAB:nj_subsetGrid:Nargin',...
      'Incorrect number of arguments');
end


switch nargout
  case 1
    data = myData;
  case 2
    data = myData;
    grd = grd;
  otherwise, error('MATLAB:nj_subsetGrid:Nargout',...
      'Incorrect number of output arguments');
end


end







