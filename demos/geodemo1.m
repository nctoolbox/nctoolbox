% NCGEODATASET GEODEMO1

url ='http://geoport.whoi.edu/thredds/dodsC/usgs/vault0/models/examples/bora_feb.nc';
nc = ncgeodataset(url)

% To access the properties we can use typical dot notation like with ordinary Matlab structures. Here we want to get a list of the variables in the
% dataset we are looking at.

 nc.variables

% The size method is a method of ncgeodataset that returns the length of each of the dimensions of a given variable in the dataset. This is a lot like
% Matlab's internal size command, but in this case we haven't even loaded any data into memory yet. All this information comes from the netcdf-java
% cdm.

 nc.size('salt')

% In this example we create a geovariable object from the salt variable in this dataset. This is done by calling geovariable with the name of the
% netcdf variable we are interested in as an argument.

 salt = nc.geovariable('salt')

% Now we can use Matlab style array indexing to subset the salt variable by its indices. Here we will assume that the arrangement of the dimensions
% follows the order of time, vertical level, horizontal coordinates. A subset of (1, 1, :, :) means that we are grabbing the first time step, the
% first level, and the entire spatial domain of the dataset.

% Note that the values that the toolbox returns are typically the same type that they are stored as in the netcdf file so they may need to be
% converted to Matlab's double.

 salinity = salt.data(1, 1, :, :);
 size(salinity)
 class(salinity)
 
% Also, it may be necessary to remove singleton dimensions for Matlab commands like plotting using the function squeeze.
 salinity= squeeze(double(salinity));

% In order to plot the salt values for the first time step/first level in the dataset as a Matlab pcolor plot, we need the spatial coordinates
% associated with the salt values. We could grab the lat and lon coordinates in the same manner that we did with the salt variable or if the data is
% CF/COARDS complaint we can take advantage of the netcdf-java common data model.

% The grid method for the geovariable object is designed to grab the all the coordinates associated with the geovariable for the given indices. Usage
% is just like the data method, except the result is a Matlab structure containing fields for each of the geovariable's dimensions and whose values
% have been subset appropriately for the requested indices.

 salinity_coords = salt.grid(1, 1, :, :)

% A higher level option is to use the grid_interop method, which returns the dimensions of our geovariable using the standardized names of lat, lon,
% time, and z instead of the original netcdf names of the coordinate dimensions.

% In addition to a more programmatic/standardized structure returned with grid_interop (interop for interoperability), the coordinate data is also
% transformed to a more standardized form:
% 
%     Coordinates recoginized by the netcdf-java cdm as time coordinates, are converted to Matlab's datenum.(This can be seen in the difference between result of the last code block and the one below.)
%     Longitude coordinates that use a 0-360 degree scheme are converted to -180 to 180 values.
%     Projected x and y values are converted to geographic coordinates lat/lon.
%     Boundary fitted vertical sigma coordinate schemes are converted to the actual vertical depth/elevation values for each grid element. (This can be seen in the difference between result of the last code block and the one below.)

     salinity_coords = salt.grid_interop(1, 1, :, :)

% Plotting using pcolor is as simple as the code below. Sometimes coordinates are stored in the netcdf datasets as vectors (vs. the 2-d arrays that
% these lat/lon coordinates are in). When this is the case, see Matlab's meshgrid function to create 2-d plaid grids from the vectors.

 pcolor(salinity_coords.lon, salinity_coords.lat, salinity)
 shading interp

% Now let's add a title to the figure that includes the dataset's global attribute title and the date of the data that we subset.

 title({nc.attribute('title'); datestr(salinity_coords.time)})