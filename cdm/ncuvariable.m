% NCUVARIABLE Extention of variable object class for unstructured grid/mesh variables.
% Trying out sort of an interface approach, I think I like it - Acrosby
%
% NCTOOLBOX (http://code.google.com/p/nctoolbox)
classdef ncuvariable < handle
    
    properties (SetAccess = private)
        dataset % matlab ncugrid class
        variable %netcdf meshvariable
        name
        axes
        size
    end
    
    methods
        
        function obj = ncuvariable(src, variableName)
            % Constructor
            if ischar(src)
                obj.dataset = ncugrid(src);  % src is a string URL/File
            elseif isa(src, 'ncugrid')
                obj.dataset = src;             % src is a ncdataset
            else
                ex = MException('NCTOOLBOX:ncuvariable', 'Invalid dataset was specified');
                ex.throw;
            end
            
            obj.variable = src.netcdfugrid.getMeshVariableByName(variableName);
            
            obj.name = variableName;
            
            coordinates = src.attribute('coordinates', variableName);
            coordinate_vars = regexp(coordinates, ' ', 'split')';
            obj.axes = coordinate_vars;
            
            obj.size = src.size(variableName);
            
        end % constructor end
        
        function s = subset(obj, struct)
            % NCUVARIABLE.SUBSET - Function to subset an unstructured model grid by lat/lon
            % bounding box using subsetting methods in the ugrid-java.
            %
            % Usage:
            %              d = uvar_obj.subset(variableName, subsetstructure)
            %
            % Subset-structure:
            %              subsetstructure.lat = [minlat maxlat];
            %              subsetstructure.lon = [minlon maxlon];
            %
            s = obj.dataset.unstructuredLatLonSubset(obj.name, struct);
        end % subset end
        
        function atts = attributes(obj)
            % NCUVARIABLE.ATTRIBUTES returns the attributes of the variable as an
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
            %     ds = ncugrid('http://somewhere/data.nc');
            %     attr = ds.attributes('time');
            %     keys = attr(:, 1);    % Cell array of keys
            %     values = attr(:, 2);  % Cell array of values
            %     i = find(ismember(keys, 'units')); % search for units attribute
            %     units = values{i};  % Retrieve the units value
            atts = obj.dataset.attributes(obj.name);
        end % attributes end
        
        function att = attribute(obj, key)
            % NCUVARIABLE.ATTRIBUTE returns the value a global attribute specified by its key or the
            % variable attribute specified by key and variable.
            %
            % Use as:
            %   a = ncdataset.attribute('title')
            %   a = ncdataset.attribute(key)
            %
            %   a = ncdataset.attribute('temp', 'title')
            %   a = ncdataset.attribute(variableName, key)
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
            att = obj.dataset.attribute(obj.name, key);
        end % attribute end
        
        function d = data(obj, first, last, stride)
            % NCUVARIABLE.DATA Retrieve all or a subset of the data for the
            % variable. The data is returned as a structure containing a
            % variable for the data as well as for each dimension of the
            % data.
            %
            % Usage:
            %   d = ncvariable.data
            %   d = ncvariable.data(first)
            %   d = ncvariable.data(first, last)
            %   d = ncvariable.data(first, last, stride)
            %
            %   If no arguments are provided all the data is returned.
            %
            % Arguments:
            %   first = The first point you want to retrieve (first point idx = 1)
            %   last  = The last point you want to retrive (default is the end of
            %       the data array)
            %   stride = The stride spacing (default is 1)
            %   NOTE! first, last, and stride must be matrices the same size as the
            %       matrix returned by NCDATASET.SIZE or SIZE
           
            d = obj.dataset.data(obj.name, first, last, stride);
        end % data end
        
        function g = grid(obj, first, last, stride)
            % NCUVARIABLE.GRID Retrieve all or a subset of the coordinate
            % data for the variable. The data is returned as a structure
            % containing a variable for each dimension of the data.
            %
            % Usage:
            %   d = ncvariable.grid
            %   d = ncvariable.grid(first)
            %   d = ncvariable.grid(first, last)
            %   d = ncvariable.grid(first, last, stride)
            %
            %   If no arguments are provided all the data is returned.
            %
            % Arguments:
            %   first = The first point you want to retrieve (first point idx = 1)
            %   last  = The last point you want to retrive (default is the end of
            %       the data array)
            %   stride = The stride spacing (default is 1)
            %   NOTE! first, last, and stride must be matrices the same size as the
            %       matrix returned by NCDATASET.SIZE or SIZE
            %
            % Returns:
            %   The data is returned as a structure containing each coordinate variable
            %

            g = obj.dataset.grid(obj.name, first, last, stride);
        end %grid end
        
        function gi = grid_interop(obj, first, last, stride)
            % NCUVARIABLE.GRID_INTEROP - Method to get the coordinate variables and their data as a
            % a structure with standardized field names for lat, lon, time, and z. Other coordiante variables
            % that are not recognized by netcdf-java as the previous four types have field names directly
            % taken from their variable names.
            % Useage: 
            %                >> gridstruct = geovar.grid_interop(1,:,:,1:2:50);
            gi = obj.dataset.grid_interop(obj.name, first, last, stride);
        end %grid_interop end
        
        %%
        function e = end(obj, k, n)
            n = obj.size;
            e = n(k);
        end % Added to deal with end indexing functionality,
        % otherwise the indexing arugment is ignored.
        
        function sref = subsref(obj,s)
            switch s(1).type
                % Use the built-in subsref for dot notation
                case '.'
                    switch s(1).subs
                        case 'grid_interop'
                            switch length(s)
                                case 1
                                    sref = obj;
                                    %                           obj.temp = size(s(2).subs{1});
                                case 2
                                    nums = obj.size;
                                    [first last stride] = indexing(s(2).subs, double(nums));
                                    sref = obj.grid_interop(first, last, stride);
                            end
                        case 'data'
                            
                            %                             if ~isempty(nums)
                            switch length(s)
                                case 1
                                    sref = obj;
                                    %                                         obj.temp = size(s(2).subs{1});
                                case 2
                                    nums = obj.size;
                                    [first last stride] = indexing(s(2).subs, double(nums));
                                    sref = obj.data(first, last, stride);
                            end
                            
                            %                             else
                            %                                 sref = obj.data;
                            %                                 warning(['NCTOOLBOX:ncgeovariable:subsref'], ...
                            %                                     ['Variable "' name '" has no netcdf dimension associated with it. Errors may result from non CF compliant files.'])
                            %                             end
                        case 'grid'
                            
                            %                             if ~isempty(nums)
                            switch length(s)
                                case 1
                                    sref = obj;
                                    %                                         obj.temp = size(s(2).subs{1});
                                case 2
                                    nums = obj.size;
                                    [first last stride] = indexing(s(2).subs, double(nums));
                                    sref = obj.grid(first, last, stride);
                            end
                            
                            %                             else
                            %                                 warning(['NCTOOLBOX:ncgeovariable:subsref'], ...
                            %                                     ['Variable "' name '" has no netcdf dimension associated with it. Errors may result from non CF compliant files.'])
                            %                             end
                        otherwise
                            sref = builtin('subsref',obj,s);
                    end
                case '()'
                    warning(['NCTOOLBOX:ncgeovariable:subsref'], ...
                        'Not a supported subscripted reference, "()" are not permitted to call variable object methods');
                case '{}'
                    warning(['NCTOOLBOX:ncgeovariable:subsref'], ...
                        'Not a supported subscripted reference, "{}" are not permitted to call variable object methods');
            end
        end % subsref end
        
    end % methods end
    
end % class end