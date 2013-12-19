%NCTOOLBOX Summary of MATLAB NCTOOLBOX capabilities.
%
%   NCTOOLBOX is a Matlab toolbox that provides read-only access to
%   common data model datasets. Under the hood, nctoolbox uses
%   NetCDF-Java as the data access layer. This allows nctoolbox to
%   access NetCDF, OPeNDAP, HDF5, GRIB, GRIB2, HDF4 and many (15+)
%   other file formats and services using the same API. It works
%   with Matlab 2008a and later.
% 
%
%   The following table lists functions supported by 
%   the nctoolbox package.
%
%      NCDATASET Functions: Basic access through NetCDF-Java
%      -----------------
%      ncdataset.size   - Return the size of a variable
%      ncdataset.data   - Fetch data
%      ncdataset.axes   - Return names of the axes of the variable
%      ncdataset.dimensions   - Returns the dimensions of the variable.
%      ncdataset.attributes   - Return the attributes of a variable
%      ncdataset.attribute    - Return an attibute of a variable
%      ncdataset.metadata     - Return the global and variable attributes.
%      ncdataset.save         - Saves the data to a local file.
%      ncdataset.time         - Converts time data to MATLAB time.
%      ncdataset.close        - Does nothing
%      ncdataset.delete       - Deletes the Return the size of a variable
%      ncvariable.size        - Return the size of a variable
%      ncvariable.attributes  - Return the attributes of a variable
%      ncvariable.attribute   - Return an attibute of a variable
%      ncvariable.axes        - Return the axes for a variable.
%      ncvariable.end         - Return the last indexes for a variable.
%      ncvariable.data        - Fetch data.
%      ncvariable.grid        - Fetch the grid associated with the data.
%
%
%      CFDATASET Functions: Access to CF-compliant extensions through NetCDF-Java
%      -----------------
%      cfdataset               - Access a CF-compliant data
%      cfdataset.standard_name - Identify a variable by standard_name
%      cfdataset.variable      - Return a ncvariable object by name
%      cfdataset.data          - Returns the data of a variable.
%      cfdataset.struct        - Returns struct data, attributes, and grid.
%      cfdataset.grid          - Return the grid of a variable
%
%      NCGEODATASET Functions: Enhanced access to CF-compliant data
%      -----------------
%      ncgeodataset               - Access a CF-compliant data
%      ncgeodataset.standard_name - Identify a variable by standard_name
%      ncgeodataset.geovariable   - Return a ncgeovariable object by name
%      ncgeodataset.data          - Returns the data of a variable.
%      ncgeodataset.struct        - Returns struct data, attributes, and grid.
%      ncgeodataset.grid          - Return the grid of a variable
%      ncgeodataset.extent        - Return the grid of a variable
%      ncgeodataset.timeextent    - Return the grid of a variable
%      ncgeodataset.gettimevar    - Return ncgeovariable of variable's time axis
%      ncgeodataset.getlonvar     - Return ncgeovariable of variable's lon axis
%      ncgeodataset.getlatvar     - Return ncgeovariable of variable's lat axis
%      ncgeodataset.gettimename   - Return name of variable's time axis
%      ncgeodataset.getlonname    - Return name of variable's lon axis
%      ncgeodataset.getlatname    - Return name of variable's lat axis
%
%      ncgeovariable              - Return a ncgeovariable object
%      ncgeovariable.axes_info    - Return a cell array of axis info
%      ncgeovariable.extent        - Return the grid of a variable
%      ncgeovariable.timeextent    - Return the grid of a variable
%      ncgeovariable.gettimevar    - Return ncgeovariable of variable's time axis
%      ncgeovariable.getlonvar     - Return ncgeovariable of variable's lon axis
%      ncgeovariable.getlatvar     - Return ncgeovariable of variable's lat axis
%      ncgeovariable.gettimename   - Return name of variable's time axis
%      ncgeovariable.getlonname    - Return name of variable's lon axis
%      ncgeovariable.getlatname    - Return name of variable's lat axis
%      ncgeovariable.getxname      - Return name of variable's x axis
%      ncgeovariable.getyname      - Return name of variable's y axis
%      ncgeovariable.geovariable   - Return a ncvariable object by name
%      ncgeovariable.data          - Returns the data of a variable.
%      ncgeovariable.grid          - Return the grid of a variable
%      ncgeovariable.grid_interop  - Return the translated grid of a variable 
%      ncgeovariable.timewindow    - Return subset of time indices
%      ncgeovariable.timewindowij  - Return the grid of a variable
%      ncgeovariable.getvectors    - Return a rotated and regridded variable 
%      ncgeovariable.geosubset     - Return a geographic subset of a variable
%      ncgeovariable.getaxesorder  - Return order of the axes variables
%
% See also netcdf, SNCTOOLS.
