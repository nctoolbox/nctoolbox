function geodemo_1a
%% NCGEODATASET GEODEMO_1A
% Method A: Read surface salinity using geovariable syntax. 
% Takes some extra steps since you create the geovariable object before 
% extracting data from it, but you then have all the geovariable 
% methods available to you.

% OPeNDAP Data URL for a not-quite CF-Compliant curvilinear ROMS model dataset
url ='http://barataria.tamu.edu:8080/thredds/dodsC/txla_nesting6/ocean_his_0196.nc';

%this data fails to be CF compliant in that its coordinate variables do not map completely to spatial coordinates

nc = ncgeodataset(url)

%% Take a look at the variables available within the dataset
% To access the properties we can use typical dot notation like with
% ordinary Matlab structures. Here we want to get a list of the variables
% in the dataset we are looking at.

 nc.variables

%% Determine the shape of the selected variable
% The size method is a method of ncgeodataset that returns the length of
% each of the dimensions of a given variable in the dataset. This is a lot
% like Matlab's internal size command, but in this case we haven't even
% loaded any data into memory yet. All this information comes from the
% netcdf-java cdm.

 nc.size('salt')

%% Create |geovariable| object to access data in MATLAB style indexing
% In this example we create a geovariable object from the salt variable in
% this dataset. This is done by calling geovariable with the name of the
% netcdf variable we are interested in as an argument.

 salt = nc.geovariable('salt')

% Now we can use Matlab style array indexing to subset the salt variable by
% its indices.  We can take a look at the dimension names using the "dimensions" 
% method:

  nc.dimensions('salt')
  
% We see the arrangement of the dimensions follows the order of time, 
% vertical level, y, and x. 
% A subset of (1, end, :, :) means that we are grabbing the first time step,
% the last vertical level, and the entire spatial domain of the dataset.

% Note that the values that the toolbox returns are typically the same type
% that they are stored as in the netcdf file so they may need to be
% converted to Matlab's double.

 salinity = salt.data(1, end, :, :);
 size(salinity)
 class(salinity)
 
% Also, it may be necessary to remove singleton dimensions for Matlab
% commands like plotting using the function squeeze.
 salinity= squeeze(double(salinity));

%% Use |geovariable| object to access coordinate data using MATLAB style indexes
% In order to plot the salt values for the first time step/first level in
% the dataset as a Matlab pcolor plot, we need the spatial coordinates
% associated with the salt values. We could grab the lat and lon
% coordinates in the same manner that we did with the salt variable or if
% the data is CF/COARDS complaint we can take advantage of the netcdf-java
% common data model.

% The grid method for the geovariable object is designed to grab the all
% the coordinates associated with the geovariable for the given indices.
% Usage is just like the data method, except the result is a Matlab
% structure containing fields for each of the geovariable's dimensions and
% whose values have been subset appropriately for the requested indices.

 salinity_coords = salt.grid(1, end, :, :)

%% Use |grid_interop| method to return coordinate axes with standard, interoperable names
% A higher level option is to use the grid_interop method, which attempts to 
% returns the dimensions of our geovariable using the standardized names of lat,
% lon, time, and z instead of the original netcdf names of the coordinate
% dimensions when possible.

% In addition to a more programmatic/standardized structure returned with
% grid_interop (interop for interoperability), the coordinate data is also
% transformed to a more standardized form:
% 
% - Time coordinates are converted to Matlab's datenum. 
% - Longitude coordinates that use a 0-360 degree scheme are converted 
%   to the range [-180, 180].
% - Projected x and y values are converted to geographic
%     coordinates lat/lon if possible.  
% - Stretched vertical coordinates are converted z values. 

 salinity_coords = salt.grid_interop(1, end, :, :)

% note that the spatial salinity coordinates are untransformed, while the time variable is
% transformed.  The variables in this file are on a staggered grid, with the bulk properties
% such as salinity or temperature on the (x_rho,y_rho,s_rho) points, and the fluxes on the (x_u,y_u) 
% or (x_v,y_v) points.  In this case, the vertical coordinate remains in non-dimensional 's' coordinates.

salinity_coords

%% Plot using MATLAB's pcolor
pcolor(salinity_coords.x_rho, salinity_coords.y_rho, salinity)
shading flat; colorbar; caxis([35 39]);

% Add a title using the global
% attribute 'title' and the date from our coordinate structure.  
 title({nc.attribute('title'); datestr(salinity_coords.time)})
