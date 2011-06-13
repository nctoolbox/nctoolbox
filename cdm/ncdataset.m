% NCDATASET  Provide access to datasets accessable by the NetCDF 4 API
%
% Use as:
%   ds = ncdataset(dataref)
%
% Arguments:
%   dataref = A reference to a ncdataset that can be accessed by the NetCDF 4
%       API. This includes local netcdf files, netcdf files on web servers
%       and OpenDAP URLs
%
% Return:
%   An instance of a ncdataset class
%
% Properties:
%   netcdf = For power users. This is an instance of
%       a ucar.nc2.ncdataset.NetcdfDataset (NetCDF-Java 4.2) and
%       is used for the underlying data access. This object can
%       be tweaked as needed. (For example, to enable/disable
%       data caching.) See
%       http://www.unidata.ucar.edu/software/netcdf-java/v4.2/javadoc/index.html
%   variables = A cell array of all variables in the ncdataset. These names
%        are used as arguments into other ncdataset methods
%
% Methods:
%   ncdataset.axes - access coordinate variable names for a given variable
%   ncdataset.attributes - access global or variable attributes
%   ncdataset.data - retrieve data (or a subset of data) for a variable
%   ncdataset.size - returns the size of a variable in the data store
%   ncdataset.time - Attempt to convert a variable to matlabs native time format (see datenum)
%
% For more information on the methods use help. For example:
%   >> help ncdataset.data
%
% Example:
%   ds = ncdataset('http://dods.mbari.org/cgi-bin/nph-nc/data/ssdsdata/deployments/m1/200810/m1_metsys_20081008_original.nc')
%   ga = ds.attributes;       % Global Attributes
%   sv = 'SonicVelocity';     % A variable that we're interested in.
%   d = ds.data(sv);          % Data for the SonicVelocity variable
%   svAx = ds.axes(sv);       % Coordinate Variable names for the SonicVelocity variable
%   svAt = ds.attributes(sv); % Attributes for SonicVelocity

% Brian Schlining (brian@mbari.org)
% 2009-05-12
% Alexander Crosby 2010, 2011
% NCTOOLBOX (http://code.google.com/p/nctoolbox)
classdef ncdataset < handle
    
    properties (SetAccess = private)
        location        % This is the location of the dataset given by the user
        netcdf          % ucar.nc2.ncdataset.NetcdfDataset java instance
        variables       % cell array containing the variable names as strings in netcdf
    end
    
    
    methods
        
        %%
        function obj = ncdataset(url)
            % NCDATASET  Constructor. Instantiates a NetcdfDataset pointing to the
            % datasource specified by 'url' and uses that as the underlying
            % dataaccess API. When instantiated, the names of all variables
            % are fetched and stored in the 'variables' property. This can be
            % use to open local files, files stored on an HTTP server and
            % OpenDAP URLs.
            %
            % If the argument is another instance of ncdataset or one of it's
            % subtypes, such as cfdataset or ncgeodataset then this creates
            % a new ncdataset that shares the underlying data.
            if isa(url, 'ncdataset')
                obj.location = url.location;
                obj.netcdf = ucar.nc2.dataset.NetcdfDataset.openDataset(url.netcdf.getLocation());
                obj.variables = url.variables;
            elseif isa(url, 'char')
                try
                    obj.netcdf = ucar.nc2.dataset.NetcdfDataset.openDataset(url);
                    
                    % Fetches the names of all variables from the netcdf ncdataset
                    % and stores them in a cell array
                    vars = obj.netcdf.getVariables();
                    n = size(vars);
                    obj.variables = cell(n, 1);
                    for i = 1:(n)
                        obj.variables{i, 1} = char(vars.get(i - 1).getName());
                    end
                    
                    obj.location = url;
                    
                catch me
                    ex = MException('NCTOOLBOX:ncdataset', ['Failed to open ' url]);
                    ex = ex.addCause(me);
                    ex.throw;
                end
            else
                ex = MException('NCTOOLBOX:ncdataset', ['The argument was a ' class(url) ...
                    '. That is not valid. Use a url (A Matlab char) or another ncdataset' ...
                    ' object instead. ']);
                ex.throw;
            end
            
            
        end
        
        
        
        %%
        function s = size(obj, variable)
            % NCDATASET.SIZE  Returns the size of the variable in the persistent store
            % without fetching the data. Helps to know what you're getting yourself
            % into. ;-)
            %
            % Use as:
            %   ds = ncdataset('http://dods.mbari.org/cgi-bin/nph-nc/data/ssdsdata/deployments/m1/200810/m1_metsys_20081008_original.nc')
            %   s = ds.size(variableName)
            %
            % Arguments:
            %   variableName = The name of the variable who's size you wish to query
            %
            % Return:
            %   The size of the data for the variable. Includes all dimensions,
            %   even the singleton dimensions
            v = obj.findvariable(variable);
            s = v.getShape()';
            
        end
        
        %%
        function d = data(obj, variable, first, last, stride)
            % NCDATASET.DATA  Fetch the data for a given variable. This methods
            % also allows you to fetch subsets of the data by specifiying a
            % the index of the first and last points as well as a stride (or
            % spacing)
            %
            % Use as:
            %   ds = ncdataset('http://dods.mbari.org/cgi-bin/nph-nc/data/ssdsdata/deployments/m1/200810/m1_metsys_20081008_original.nc')
            %   d = ds.data(variableName)
            %   d = ds.data(variableName, first)
            %   d = ds.data(variableName, first, last)
            %   d = ds.data(variableName, first, last, stride)
            %
            % Arguments:
            %   variableName = The name of the variable whos data you want to retrieve
            %   first = The first point you want to retrieve (first point idx = 1)
            %   last  = THe last point you want to retrive (default is the end of
            %       the data array)
            %   stride = The stride spacing (default is 1)
            %   NOTE! first, last, and stride must be matrices the same size as the
            %       matrix returned by NCDATASET.SIZE
            %
            % Return:
            %   The data from the variable. The data will be returned in the
            %   matlab type that nearest corresponds to the variables data
            %   type in the data store. For most Matlab operations you will
            %   need to convert the data to double precision. Here's an example:
            %     ds = ncdataset('/Volumes/oasis/m1/netcdf/m1_adcp_20051020_longranger.nc')
            %     u = ds.data('u_component_uncorrected'); % u is a matlab 'single' type
            %     u = double(u) % promote single to double precision
            switch nargin
                case 2
                    d = obj.readdata(variable);
                case 3
                    d = obj.readdata(variable, first);
                case 4
                    d = obj.readdata(variable, first, last);
                case 5
                    d = obj.readdata(variable, first, last, stride);
            end
            
        end
        
        
        %%
        function cv = axes(obj, variable)
            % NCDATASET.AXES  Returns a cell array containing the variable names of
            % coordinate axes for the given variable
            %
            % Use as:
            %   ds = ncdataset('http://dods.mbari.org/cgi-bin/nph-nc/data/ssdsdata/deployments/m1/200810/m1_metsys_20081008_original.nc')
            %   ax = ds.axes(variableName);
            %
            % Inputs:
            %   variableName = the name of the variable whos axes you want to retrieve
            %
            % Return:
            %   An (n, 1) cell array containing the names (in order) of the variable
            %   names representing the coordinate axes for the specified variableName.
            v = obj.findvariable(variable);
            dims = v.getDimensions();
            cv = cell(dims.size(), 1);
            for i = 1:dims.size();
                ca = obj.netcdf.findCoordinateAxis(dims.get(i - 1).getName());
                if ~isempty(ca)
                    cv{i} = char(ca.getName());
                end
            end
        end
        
        %%
        function a = attributes(obj, variable)
            % NCDATASET.ATTRIBUTES returns the attributes of the variable as an
            % n x 2 cell array.
            %
            % Use as:
            %   a = ncdataset.attributes
            %   a = ncdataset.attributes(variableName)
            %
            % Inputs:
            %   variableName = The name of the variable whos attributes you want
            %       to retrieve. If no argument is specified then the
            %       global attributes are returned.
            %
            % Return:
            %   An (n, 2) cell array. The 1st column contains the attribute name
            %   The 2nd column contains the attribute value. Each row is essentially
            %   a key-value pair.
            %
            % Hints:
            %   Here's how to obtain cell arrays of the keys and corresponding values
            %   and search for a particular key (or attribute name)
            %     ds = ncdataset('http://somewhere/data.nc');
            %     attr = ds.attributes('time');
            %     keys = attr(:, 1);    % Cell array of keys
            %     values = attr(:, 2);  % Cell array of values
            %     i = find(ismember(keys, 'units')); % search for units attribute
            %     units = values{i};  % Retrieve the units value
            if (nargin < 2)
                % Get global attributes
                aa = obj.netcdf.getGlobalAttributes();
            else
                v = obj.findvariable(variable);
                aa = v.getAttributes();
            end
            
            if ~isempty(aa)
                n = aa.size();
                a = cell(n, 2);
                for i = 1:n
                    at = aa.get(i - 1);
                    a{i, 1} = char(at.getName());
                    if (at.isString())
                        a{i, 2} = char(at.getStringValue());
                    else
                        a{i, 2} = at.getValues().copyToNDJavaArray();
                    end
                end
            else
                % Show warning, return empty cell array
                warning('NCTOOLBOX:ncdataset:attributes', 'No attributes were found');
                a = cell(1, 2);
            end
        end
        
        function val = attribute(obj, key, var)
            % NCDATASET.ATTRIBUTE returns the value a global attribute specified by its key or the
            % variable attribute specified by key and variable.
            %
            % Use as:
            %   a = ncdataset.attribute('title')
            %   a = ncdataset.attribute(key)
            %
            %   a = ncdataset.attribute('title', 'temp')
            %   a = ncdataset.attribute(key, variableName)
            %
            % Inputs:
            %   key = The name of the attribute field like 'title' or 'units'...
            %   variableName = The name of the variable whos attributes you want
            %       to retrieve. If no argument is specified then the
            %       global attributes are returned.
            %
            % Return:
            %   The value associated with the attribute field corresponding to key (and optionally
            %   variableName)
            if nargin < 3
                atlist = obj.attributes;
                val = value4key(atlist, key);
            elseif nargin > 1
                atlist = obj.attributes(var);
                val = value4key(atlist, key);
            else
                warning('NCTOOLBOX:ncdataset:attribute', 'No key or variable specified.');
            end
        end
        %%
        function save(obj, filename)
            % NCDATASET.SAVE Save the data to a local netcdf file
            %
            % Use as:
            %   obj.save(filename)
            %
            % Arguments:
            %   filename = The path, as a string, to the file that the data will
            %       be written to.
            
            % TODO Should we warn user if file already exists?
            ucar.nc2.FileWriter.writeToFile(obj.netcdf, filename);
        end
        
        %%
        function t = time(obj, variable, data)
            % NCDATASET.TIME  Attempts to convert data to Matlab's native time format
            %
            % Use as:
            %   t = obj.time(variableName)
            %   t = obj.time(variableName, data)
            %
            % Arguments:
            %   variableName = The name of the variable you're trying to convert.
            %       This is needed in order to fetch information about the variable.
            %       Currently, this method depends on the presence of a
            %       variable attribute named 'units'.
            %   data = the data you're trying to convert. If not provided then
            %       all the data for the given variableName is fetched and converted.
            %       It's useful to supply the data when you are working with
            %       a subset of the data.
            
            if nargin < 3
                data = obj.data(variable);
            end
            try
                t = convertToTime1(obj, variable, data);
            catch me
                try
                    t = convertToTime2(obj, variable, data);
                catch me2
                    warning('NCTOOLBOX:ncdataset:time', ['Unable to convert ' variable ' to matlab time.']);
                    t = data;
                end
            end
            
        end
        
        %%
        function close(obj)
            % NCDATASET.CLOSE Added to prevent throwing an error. This is problematic since close is a 
            % keyword in matlab for closing figure windows, but requested by Rich to avoid throwing errors
            % in legacy njtbx code...
        end
        
        function delete(obj)
            % NCDATASET.DELETE Closes netcdf files when object NCDATASET object is disposed or leaves scope
            try
                obj.netcdf.close()
            catch me
                % Do nothing
            end
        end
        
        % TODO Add ncml method that returns a string of NcML for this
        % netcdf dataset
        
    end
    methods (Access = protected)
        
        function v = findvariable(obj, variable)
            % NETCDF.FINDVARIABLE - Helper function that will escape a variable name if needed. 
            v = obj.netcdf.findVariable(variable);
            if isempty(v)
                v = obj.netcdf.findVariable(ucar.nc2.NetcdfFile.escapeName(variable))
            end
            if isempty(v)
                warning('NCTOOLBOX:ncdataset:findvariable', ['Could not find the variable: ' variable]);
            end
        end
        
        %%
        function d = readdata(obj, variable, first, last, stride)
            % NETCDF.READDATA - Helper function that's called by NETCDF.DATA
            v = obj.findvariable(variable);
            
            if (nargin == 2)
                array = v.read();
                try
                    d = array.copyToNDJavaArray(); % this fails if the variable has no java shape/no dimension was assigned
                catch me1
                    try
                        % TODO (Alex added this code) Where is a file where
                        % this code section gets called?
                        d = array.toString;  % different way to get single value out of java array
                        d = d.toCharArray';  % must transpose
                        d = str2double(d);   % matlab string to matlab numeric
                    catch me2
                        ex = MException('NCTOOLBOX:ncdataset:data', ['Failed to open "' variable '" in ' url]);
                        ex = ex.addCause(me2);
                        ex.throw;
                    end
                end
                %           d = double(d);
            else
                s = obj.size(variable);
                
                % Fill in missing arguments
                % default stride is 1
                if (nargin < 5)
                    stride = ones(1, length(s));
                end
                
                % Default last is the end
                if (nargin < 4)
                    last = s;
                end
                
                % Construct the range objects needed to subset the data
                n = max(size(obj.size(variable)));
                ranges = java.util.ArrayList(n);
                for i = 1:n
                    ranges.add(ucar.ma2.Range(first(i) - 1, last(i) - 1, stride(i)));
                end
                
                array = v.read(ranges);
                d = array.copyToNDJavaArray();
                
            end
            
        end
    end
    
end

%%
function d = convertToTime1(obj, variable, data)
% convertToTime1 - Attempt time units conversion using udunits
attr = obj.attributes(variable);
units = value4key(attr, 'units');
dateUnit = ucar.nc2.units.DateUnit(['0 ' units]);
[r c] = size(data);
d = ones(r, c) * NaN;

% Can't vectorize this operation, use slow-ass loop instead.
for i = 1:r
    for j = 1:c
        d(i, j) = utc2sdn(dateUnit.makeDate(data(i, j)).getTime() / 1000);
    end
end

end

%%
function t = convertToTime2(obj, variable, data)
% convertToTime2 - Attempt time units conversion by hand

% NOTE: THis is naive implementation that looks for  aunit attribute on the
% variable with a value in the form of '[time unit] since [some date]'
attr = obj.attributes(variable);
keys = attr(:, 1);
values = attr(:, 2);
i = find(ismember(keys, 'units'));         % search for units attribute
units = lower(values{i});                  % Retrieve the units value


% Figure out the conversion from the time units to days (matlabs time units)
conversion = 1;
unitPatterns = {'milli', 'sec', 'minute', 'hour', 'day', 'year'};
conversions = [1 / (1000 * 60 * 60 * 24), 1 / (60 * 60 * 24), 1 / (60 * 24), 1 /24, 1, 365.25];
for i = 1:length(unitPatterns)
    timeUnit = unitPatterns{i};
    u = strfind(units, timeUnit);
    if ~isempty(u)
        conversion = conversions(i);
        break
    end
end

% Figure out the offset
offset = 0;
idx = regexp(units, '\d');
if ~isempty(idx)
    startIdx = idx(1);
    endIdx = idx(end);
    dateString = units(startIdx:endIdx);
    % Try using matlabs date parsing which covers most common cases
    try
        offset = datenum(dateString);
    catch e
        % If we're here then it's most likey an iso8601 format. We'll try one.
        try
            df = java.text.SimpleDateFormat('yyyy-MM-dd''t''HH:mm:ss');
            df.setTimeZone(java.util.TimeZone.getTimeZone('UTC'));
            jDate = df.parse(dateString);
            secs = jDate.getTime() / 1000;
            offset = utc2sdn(secs);
        catch me
            warning('NCTOOLBOX:ncdataset:convertToTime2', ['Unable to parse date: ' dateString]);
        end
    end
end


if (conversion == 1) && (offset == 0)
    warning('NCTOOLBOX:ncdataset:convertToTime2', ['No conversion occurred. Are you sure that ' variable ' is time data?']);
end

% fprintf(1, 'Conversion = %12.5f, Offset = %12.5f\n', conversion, offset);
t = data * conversion + offset;

end




