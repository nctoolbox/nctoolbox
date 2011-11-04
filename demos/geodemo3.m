% NCGEODATASET GEODEMO2

url ='http://geoport.whoi.edu/thredds/dodsC/usgs/vault0/models/examples/bora_feb.nc';
nc = ncgeodataset(url)

% To access the properties we can use typical dot notation like with ordinary Matlab structures. Here we want to get a list of the variables in the
% dataset we are looking at.

 nc.variables

% The size method is a method of ncgeodataset that returns the length of each of the dimensions of a given variable in the dataset. This is a lot like
% Matlab's internal size command, but in this case we haven't even loaded any data into memory yet. All this information comes from the netcdf-java
% cdm.

 nc.size('salt')

% This syntax should be familiar to those who have seen njTBX or mexcdf toolboxes for Matlab. In fact, the function for this syntax is intended to
% explicitly replicate the functionality of the same syntax in njTBX to help transition users and their legacy code.

% To return the values of the salt variable simply add the variable name to the braces argument, and the Matlab style indices subset to the
% parenthesis argument. (If no indices are included, a geovariable object is returned. (:) can be used for the entire variable.)

 salinity = nc{'salt'}(1, 1, :, :);
 size(salinity)

 salinity = squeeze(double(salinity));

% The grid method in this syntax is the exact same command as the grid_interop method in the approach from Method 1. As in accessing the variable
% values, use the same arguments to access the coordinate values for the given variable and subset, but append a .grid at the end.

% The coordinate data is also transformed to a more standardized form:
% 
%     Coordinates recoginized by the netcdf-java cdm as time coordinates, are converted to Matlab's datenum.
%     Longitude coordinates that use a 0-360 degree scheme are converted to -180 to 180 values.
%     Projected x and y values are converted to geographic coordinates lat/lon.
%     Boundary fitted vertical sigma coordinate schemes are converted to the actual vertical depth/elevation values for each grid element. 

 salinity_coords = nc{'salt'}(1, 1, :, :).grid

% Plotting using pcolor is as simple as the code below. Sometimes coordinates are stored in the netcdf datasets as vectors (vs. the 2-d arrays that
% these lat/lon coordinates are in). When this is the case, see Matlab's meshgrid function to create 2-d plaid grids from the vectors.

 pcolor(salinity_coords.lon, salinity_coords.lat, salinity)
 shading interp

% Now let's add a title to the figure that includes the dataset's global attribute title and the date of the data that we subset.

 title({nc.attribute('title'); datestr(salinity_coords.time)})