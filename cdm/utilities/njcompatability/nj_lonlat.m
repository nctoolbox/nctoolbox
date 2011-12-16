function [lon,lat]=nj_lonlat(ncRef,var);
% NJ_LONLAT  find longitude and latitude variables
% [lon,lat]=nj_lonlat(fname,var);
%    Inputs: fname = file name, url, or netcdf file object
%            var   = variable name (e.g. 'U10')
%         idouble  = [0] default, [1] convert output to double
%    Output: lon   = longitude array corresponding to "var"
%            lat   = latitude array corresponding to "var"

% import the NetCDF-Java methods we need for this example

% Rich Signell (rsignell@usgs.gov);

if nargin < 2, help(mfilename), return, end

try
    if (isa(ncRef, 'ncgeodataset')) %check for ncgeodataset Object
        nc = ncRef;
    else
        % open CF-compliant NetCDF File as a Common Data Model (CDM) "Grid Dataset"
        nc = ncgeodataset(ncRef);
    end
    gvar=nc.geovariable(var);
    lon=nc{gvar.getlonname}(:);
    lat=nc{gvar.getlatname}(:);
catch
    %gets the last error generated
    err = lasterror();
    disp(err.message);
end
