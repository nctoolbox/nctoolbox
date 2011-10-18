function [data,grd]=nj_tslice(ncRef,var,iTime, level)
% NJ_TSLICE - Get data and coordinates from a CF-compliant file at a specific time step and level
%
% Usage:
%   [data,grd]=nj_tslice(ncRef,var,[iTime], [level]);
% where,
%   ncRef - Reference to netcdf file. It can be either of two
%           a. local file name or a URL  or
%           b. An 'ncgeodataset' matlab object, which is the reference to already
%              open netcdf file.
%              [ncRef=ncgeodataset(uri)]
%   var - variable to slice
%   iTime - time step to extract data. Use 'inf' for last time step.
%   level - vertical level (if not supplied, volume data is retrieved.) Use
%   'inf' for last level.
% returns,
%   data - data  - matlab array
%          grd   - structure containing lon,lat,z,jdmat (Matlab datenum)
% e.g.,
%   uri ='http://geoport.whoi.edu/thredds/dodsC/usgs/vault0/models/examples/bora_feb.nc';% NetCDF/OpenDAP/NcML file
%   [data,grd]=nj_tslice(uri,'temp',2, 14); % Retrieve data from level 14 at time step 2
%   [data,grd]=nj_tslice(uri,'temp',2); % Retrieve 3D data at time step 2
%   [data,grd]=nj_tslice(uri,'h'); % Retrieve 2D non time dependent array
%
%
% rsignell@usgs.gov
% Sachin Kumar Bhate (skbhate@ngi.msstate.edu)

if nargin < 2, help(mfilename), return, end

%initialize 
% data=[];
% grd.lat=[];
% grd.lon=[];
% grd.z=[];
% grd.time=[];
isNcRef=0;

try
    if (isa(ncRef, 'ncgeodataset')) %check for ncgeodataset Object
        nc = ncRef;
        isNcRef=1;
    else
        % open CF-compliant NetCDF File as a Common Data Model (CDM) "Grid Dataset"
        nc = ncgeodataset(ncRef);
    end
    
    gvar = nc.geovariable(var); %get geovariable
    
    switch nargin
        case 2
            % read the volume data (3D). All times.
            data = squeeze(gvar.data(:,:,:,:));
            if nargout==2,
                grd = gvar.grid_interop(:,:,:,:);
            end
        case 3
            if (isinf(iTime) || iTime == -1)
                data = squeeze(gvar.data(end,:,:,:));
                if nargout==2,
                    grd = gvar.grid_interop(end,:,:,:);
                end
            else
                data = squeeze(gvar.data(iTime,:,:,:));
                if nargout==2,
                    grd = gvar.grid_interop(iTime,:,:,:);
                end
            end
        case 4
            if ( (isinf(iTime)) || (isinf(level)) )
                a=size(tvar);
                if isinf(level);level=a(2);end
                if isinf(iTime);iTime=a(1);end
            end
            data = squeeze(gvar.data(iTime,level,:,:));
            if nargout==2,
                grd = gvar.grid_interop(iTime,level,:,:);
            end

    end

catch
    %gets the last error generated
    err = lasterror();
    disp(err.message);
end

