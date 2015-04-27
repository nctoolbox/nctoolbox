function setnccache(cachepath)
% SETNCCACHE Configure the NetCDF cache location.
%
% Usage:
%    setnccache(cachepath)
%
% Inputs:
%   cachepath: The directory for NetCDF to write temporary files into, such as
%       grib index files.

% Brian Schlining
% 2015-04-27

import ucar.nc2.grib.GribCollection;
import ucar.nc2.util.DiskCache;
import ucar.nc2.util.DiskCache2;

% Inform NJ where the cache is and to always use it.
DiskCache.setRootDirectory(cachepath);
DiskCache.setCachePolicy(true);

% The following is necessary as of NJ 4.3.15 in order to prevent GRIB
%temp files from cluttering up same directory as source GRIB files.
dc2 = DiskCache2(cachepath, false, 30, 30);
dc2.setAlwaysUseCache (true);
GribCollection.setDiskCache2 (dc2);

end
