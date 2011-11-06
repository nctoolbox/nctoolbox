function [data,grd]=nj_tslice(ncRef,var,iTime,level)
% NJ_TSLICE - Get data and coordinates from CF-compliant dataset at specific time and level
% Usage:
%   [data,grd]=nj_tslice(ncRef,var,[iTime],[level]);
% Inputs: 
%   ncRef = OpenDAP Data URL, Local file, or existing 'ncgeodataset' object
%   var - variable to slice
%   iTime = time step to extract data. Use inf or -1 for last time step.
%   level = vertical level (if not supplied, volume data is retrieved). Use 
%           inf or -1 for last level.
% Outputs:
%   data = data  - matlab array
%   grd  = structure containing lon,lat,z,time (Matlab datenum)
% Example:
%   uri ='http://geoport.whoi.edu/thredds/dodsC/examples/bora_feb.nc';% NetCDF/OpenDAP/NcML file
%   [data,grd]=nj_tslice(uri,'temp',2, 14); % get level 14 at time step 2
%   [data,grd]=nj_tslice(uri,'temp',2); % get 3D data at time step 2
%   [data,grd]=nj_tslice(uri,'h'); % Retrieve 2D non time dependent array
%
% NCTOOLBOX (http://code.google.com/p/nctoolbox)

if nargin < 2, help(mfilename), return, end
try
  if (isa(ncRef, 'ncgeodataset')) %check for ncgeodataset object
    nc = ncRef;
  else
    % open CF-compliant dataset as ncgeodataset object
    nc = ncgeodataset(ncRef);
  end 
  gvar = nc.geovariable(var); %get geovariable object
  switch nargin
    case 2  % return 2D or 3D data (z,lat,lon) for all times (could be huge)
      data = squeeze(gvar.data(:,:,:,:));
      if nargout==2,
        grd = gvar.grid_interop(:,:,:,:);
        grd.z = squeeze(grd.z);
      end
    case 3  % return 2D or 3D data at a specific time 
      a=size(gvar);
      if (isinf(iTime) || iTime == -1)
        iTime=a(1);  % assumes time is 1st dimension
      end
      data = squeeze(gvar.data(iTime,:,:,:));
      if nargout==2,
        grd = gvar.grid_interop(iTime,:,:,:);
        grd.z = squeeze(grd.z);
      end
    case 4 % return a single 2D at a specific time and layer  
      a=size(gvar);
      if (isinf(iTime) || iTime == -1)
        iTime=a(1);  % assumes time is 1st dimension
      end
      if (isinf(level) || level == -1)
        level=a(2);  % assumes level is 2nd dimension 
      end
      data = squeeze(gvar.data(iTime,level,:,:));
      if nargout==2,
        grd = gvar.grid_interop(iTime,level,:,:);
        grd.z = squeeze(grd.z);
      end     
  end 
catch ME
  ME.message
end

