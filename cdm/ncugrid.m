% NCUGRID Extention of dataset object class for unstructured grid datasets.
% Usage:
%            >> dataset = ncugrid(uri);
%
% NCTOOLBOX (http://code.google.com/p/nctoolbox)
classdef ncugrid < handle
    
    properties (SetAccess = private)
        dataset % ncdataset obj (may not even be needed) since I am reproducing essential functions in
        % this ncugrid class...
        netcdfugrid % Kyle's netcdf u grid dataset
        numberofgrids % total number of meshes
        meshes % list of mesh objs
        cells % number of cells
        nodes % number of nodes
        edges % number of edges
        faces % number of faces
        variables % list of variable names
    end
    
    methods
        
        function obj = ncugrid(nc)
            import java.io.IOException;
            import java.util.Formatter;
            import java.util.List;
            import ucar.nc2.Variable;
            import ucar.nc2.constants.FeatureType;
            import ucar.nc2.dt.UGridDataset;
            import ucar.nc2.dt.ugrid.Cell;
            import ucar.nc2.dt.ugrid.Entity;
            import ucar.nc2.dt.ugrid.Mesh;
            import ucar.nc2.dt.ugrid.geom.LatLonPoint2D;
            import ucar.nc2.ft.FeatureDatasetFactoryManager;
            import ucar.nc2.util.CancelTask;
            import ucar.unidata.geoloc.LatLonPoint;
            import ucar.unidata.geoloc.LatLonPointImpl;
            import ucar.unidata.geoloc.LatLonRect;
            
            if ischar(nc)
                try
                    obj.dataset = ncdataset(nc);  % src is a string URL/File
                catch % if it doesnt work maybe try dods instead of http
                    nc(1:4) = 'dods';
                    obj.dataset = ncdataset(nc);  % src is a string URL/File
                end
                cancelTask = [];
                obj.netcdfugrid = FeatureDatasetFactoryManager.open(FeatureType.UGRID, nc, cancelTask, Formatter());
                
            elseif isa(nc, 'ncdataset')
                obj.dataset = nc;             % src is a ncdataset
                cancelTask = [];
                obj.netcdfugrid = FeatureDatasetFactoryManager.open(FeatureType.UGRID, nc.link, cancelTask, Formatter());
            elseif isa(nc,'ucar.nc2.dt.ugrid.UGridDataset')
                obj.dataset = nc;             % src is a raw ugriddataset in memory
                cancelTask = [];
                obj.netcdfugrid = FeatureDatasetFactoryManager.open(FeatureType.UGRID, nc, cancelTask, Formatter());
            else
                ex = MException('NCUGRID:ncugrid', 'Invalid dataset was specified');
                ex.throw;
            end
            
            %set properties
            u = obj.netcdfugrid;
            
            vars = u.getMeshVariables();
                    n = size(vars);
                    obj.variables = cell(n, 1);
                    for i = 1:(n)
                        obj.variables{i, 1} = char(vars.get(i - 1).getName());
                    end
                    
            meshsets = u.getMeshsets();
            %       UGridDataset.getMeshsets().get(0).getMesh()
            obj.meshes = meshsets;
            obj.numberofgrids = meshsets.size;
            for i = 1:(obj.numberofgrids)
                obj.cells(i, 1) = meshsets.get(i - 1).getMesh.getSize();
                obj.nodes(i, 1) = meshsets.get(i - 1).getMesh.getNodeSize();
                obj.edges(i, 1) = meshsets.get(i - 1).getMesh.getEdgeSize();
                obj.faces(i, 1) = meshsets.get(i - 1).getMesh.getFaceSize();
            end
            
            
        end % end constructor
        
        function uvar = uvariable(obj, varName)
            % NCUGRID.UVARIABLE - Instantiates a ugrid variable object from the variable 'varName'.
            % Usage:
            %             uvar_obj = nc.uvariable(varName);
            %
            uvar = ncuvariable(obj, varName);
        end % uvariable end
        
        function ss = unstructuredLatLonSubset(obj, varargin) % input subset structure
            % NCUGRID.unstructuredLatLonSubset  - Function to subset an unstructured model grid by lat/lon
            % bounding box using subsetting methods in the ugrid-java.
            %
            % Usage:
            %              d = nc.unstructuredLatLonSubset(variableName, subsetstructure)
            %
            % Subset-structure:
            %              subsetstructure.lat = [minlat maxlat];
            %              subsetstructure.lon = [minlon maxlon];
            %
            import ucar.nc2.dt.ugrid.geom.LatLonPolygon2D;
            import ucar.nc2.dt.ugrid.geom.LatLonRectangle2D;
            import ucar.unidata.geoloc.LatLonPoint;
            import ucar.unidata.geoloc.LatLonPointImpl;
            import ucar.unidata.geoloc.LatLonRect;
            if nargin ==3
                struct = varargin{2};
                varName = varargin{1};
                % Find the mesh that contains the variable of interest.
                mesh_ident = obj.attribute('mesh', varName);
                mesh_vars = regexp(mesh_ident, ' ', 'split');
                mesh_ind = str2double(mesh_vars{1}(5:end));
            elseif nargin == 2
                struct = varargin{1};
            end
            % Need to construct lat/lon rect:
            maxlat  = max(struct.lat);
            minlon = min(struct.lon);
            minlat = min(struct.lat);
            maxlon = max(struct.lon);
            bbox = LatLonRect(LatLonPointImpl(maxlat, minlon), LatLonPointImpl(minlat, maxlon));
            
            
            if nargin == 3
                % Get all the meshes, access the correct mesh for the variable, then call subset method on mesh.
                meshes = obj.meshes;
                mesh = meshes.get(mesh_ind-1).getMesh();
                if mesh.getTreeSize() == 0
                    mesh.buildRTree();
                end
%                 v = obj.netcdfugrid.getMeshVariableByName(varName);
%                 subsat = v.subsetToDataset(bbox); % Subsat is ugriddataset
                v = obj.netcdfugrid;
                subsat = v.subset(bbox); % Subsat is ugriddataset
            elseif nargin == 2
                v = obj.netcdfugrid;
                subsat = v.subset(bbox); % Subsat is ugriddataset
            end
           
            % Get netcdfj variable class and read all (since we already subset)
            submeshv = subsat.getMeshVariableByName(varName);
            connect_meshset = submeshv.getMeshset();
            temp = connect_meshset.getMesh();
            temp.toString()
            
            % Get variable values from subsat meshVariable
            submeshv = submeshv.getVariable();
            array = submeshv.read();
            values = array.copyToNDJavaArray();
            ss.data = values;
            % Related coordinates and connectivity
            ss.grid = obj.getSomeCoordinateData(varName, connect_meshset);
%         end
            
        end % end subset
        
        function a = attributes(obj, variable)
            % NCUGRID.ATTRIBUTES returns the attributes of the variable as an
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
            if (nargin < 2)
                % Get global attributes
                aa = obj.netcdfugrid.getGlobalAttributes();
            else
                % Get attributes for the variable
                variable = ucar.nc2.NetcdfFile.escapeName(variable);
                v = obj.netcdfugrid.getNetcdfFile().findVariable(variable);
                if isempty(v)
                    warning('NCTOOLBOX:ncugrid:attributes', ['Could not find the variable ']);
                end
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
                warning('NCTOOLBOX:ncugrid:attributes', 'No attributes were found');
                a = cell(1, 2);
            end
        end % end attributes
        
        function val = attribute(obj, varargin)
            % NCUGRID.ATTRIBUTE returns the value a global attribute specified by its key or the
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
            if nargin < 3
                atlist = obj.attributes;
                val = value4key(atlist, varargin{1});
            elseif nargin > 1
                atlist = obj.attributes(varargin{1});
                val = value4key(atlist, varargin{2});
            else
                warning('NCTOOLBOX:ncugrid:attribute', 'No key or variable specified.');
            end
        end % end attribute
        
        function gs = getSomeCoordinateData(obj, varName, meshset)
            mesh_ident = obj.attribute('mesh', varName);
            mesh_vars = regexp(mesh_ident, ' ', 'split');
            if nargin < 3
                meshset  = obj.meshes.get(str2double(mesh_vars{1}(5:end)-1));
            end
            
            coordinates = obj.attribute('coordinates', varName);
            coordinate_vars = regexp(coordinates, ' ', 'split');
            
                  location = obj.attribute('location', varName);
                  connectivity_var = obj.attribute([location, '_connectivity'], mesh_vars{1});
            
            % Get coordinates (time, lat, lon,...)
            for i = 1:length(coordinate_vars)
                switch coordinate_vars{i}
                    case 'time'
                        v = meshset.getMeshVariableByName(coordinate_vars{i});
                        if isempty(v)
                            v = obj.netcdfugrid.getNetcdfFile().findVariable(coordinate_vars{i});
                        else
                            %
                        end
                        array = v.read();
                        values = array.copyToNDJavaArray();
                        gs.(coordinate_vars{i}) = values;
                        clear values
                    otherwise
                        
                end
            end % end for
            
            nodes = meshset.getMesh().getUniqueNodes();
            % Kyle made a method to do this automatically mesh.getNodeLatLons() and getNodeIndexes()
            for i = 1:nodes.size();
                point = nodes.get(i-1).getGeoPoint();
                gs.lat(i,1) = point.getLatitude();
                gs.lon(i,1) = point.getLongitude();
                gs.index(i,1) = nodes.get(i-1).getDataIndex();
            end
            
            try
                % Get connectivity
                %                     v = obj.netcdfugrid.getNetcdfFile().findVariable(connectivity_var);
                v = meshset.getMesh();
                v = v.getConnectivityVariable();
                v = v.getVariable();
                
                array = v.read();
                values = array.copyToNDJavaArray();
                gs.connectivity = values;
            catch
            end
            
        end % end
        
        function d = data(obj, variable, first, last, stride)
            % NCUGRID.DATA - Get variable data based on vector of elements (nodes, faces, etc.) representation of the unstructured data.
            %
            % Usage:
            %              d = nc.data(variableName)
            %              d = nc.data(variableName, first)
            %              d = nc.data(variableName, first, last)
            %              d = nc.data(variableName, first, last, stride)
            %
            
            variable = ucar.nc2.NetcdfFile.escapeName(variable);
            v = obj.netcdfugrid.getNetcdfFile().findVariable(variable);
            
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
                        ex = MException('NCTOOLBOX:ncugrid:data', ['Failed to open "' variable '" in ' url]);
                        ex = ex.addCause(me2);
                        ex.throw;
                    end
                end
                d = d;
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
        end % end data
        
        function ig = grid_interop(src, varName, first, last, stride)
            % NCUGRID.GRID_INTEROP - Method to get the coordinate variables and their data as a
            % a structure with standardized field names for lat, lon, time, and z. Other coordiante variables
            % that are not recognized by netcdf-java as the previous four types have field names directly
            % taken from their variable names.
            % Usage: 
            %          >> gridstruct = geovar.grid_interop(1,:,:,1:2:50);
            %
            g = src.grid(varName, first, last, stride);
            names = fieldnames(g);
            
            for i = 1:length(names); % loop through fields returned by grid
                tempname = names{i};
                javaaxisvar  =   src.dataset.netcdf.findVariable(tempname);
                try
                    type = char(javaaxisvar.getAxisType);
                    if isempty(type)
                        ig.(tempname) = g.(tempname);
                    else
                        switch type
                            case 'Height'
                                pos_z = char(javaaxisvar.getPositive());
                                if strcmp(pos_z, 'POSITIVE_DOWN')
                                    tmp = g.(tempname);
                                    ig.z = tmp.*-1; %adjust for positive direction
                                else
                                    ig.z = g.(tempname);
                                end
                                
                            case 'GeoZ'
                                pos_z = char(javaaxisvar.getPositive());
                                if strcmp(pos_z, 'POSITIVE_DOWN')
                                    tmp = g.(tempname);
                                    ig.z = tmp.*-1; %adjust for positive direction
                                else
                                    ig.z = g.(tempname);
                                end
                                
                            case 'Time'
                                tmp = g.(tempname);
                                t_converted = src.dataset.time(tempname, tmp);
                                ig.time = t_converted;
                                
                                %                     case 'RunTime'
                                %                       tmp = obj.dataset.data(name, vFirst, vLast, vStride);
                                %                       t_converted = obj.dataset.time(name, tmp);
                                %                       data.(type) = t_converted;
                                
                            case 'Lon'
                                tmp = g.(tempname);
                                ind = find(tmp > 180); % convert 0-360 convention to -180-180
                                tmp(ind) = tmp(ind)-360;
                                ig.lon = tmp;
                                
                            case 'Lat'
                                ig.lat = g.(tempname);
                                
                            otherwise
                                ig.(tempname) = g.(tempname);
                                
                        end % end switch on type
                    end % end is type empty or not if statement
                catch
                    ig.(tempname) = g.(tempname);
                end
            end % end loop through field names
            
        end % grid_interop end
        
        function gr = grid(obj, varName, first, last, stride)
            % NCUGRID.GRID - Return a matlab structure of the coordinate variables of 'variableName',
            % subsetted using first, last, and stride index arguments in relation to the variable of interest
            % 'variableName'.
            %
            % Usage:
            %              d = nc.grid(variableName)
            %              d = nc.grid(variableName, first)
            %              d = nc.grid(variableName, first, last)
            %              d = nc.grid(variableName, first, last, stride)
            %
            mesh_ident = obj.attribute('mesh', varName);
            mesh_vars = regexp(mesh_ident, ' ', 'split');
            if nargin < 3
                mesh  = obj.meshes.get(str2double(mesh_vars{1}(5:end)-1));
            end
            
            coordinates = obj.attribute('coordinates', varName);
            coordinate_vars = regexp(coordinates, ' ', 'split');
            
            location = obj.attribute('location', varName);
            connectivity_var = obj.attribute([location, '_connectivity'], mesh_vars{1});
            
            % Get coordinates (time, lat, lon,...)
            for i = 1:length(coordinate_vars)
                
                
                v = obj.netcdfugrid.getNetcdfFile().findVariable(coordinate_vars{i});
                vs = obj.size(v.getName);
                s = obj.size(varName);
                %         if numel(vs) > 0 % Added to solve work around somedata calls that involve variables with
                %           % no netcdf dim. (This will be frequent in some OOI files.)
                if (length(vs) == length(s))
                    %% case: sizes are the same
                    if isequal(vs, s)
                        vFirst = first;
                        vLast = last;
                        vStride = stride;
                    else
                        me = MException('NCTOOLBOX:ncvariable:somedata', ...
                            ['The data size of the coordinate variable,' ...
                            name ', does not fit the size of ' obj.name]);
                        me.throw;
                    end
                    
                elseif length(vs) == 1
                    %% case: singleton dimension. Find side of data with
                    % the same length
                    
                    % TODO: the following line  will give bogus results if
                    % the data has 2 dimensions of the same length
                    dim = find(s == vs, 1);
                    if ~isempty(dim)
                        vFirst = first(dim);
                        vLast = last(dim);
                        vStride = stride(dim);
                    else
                        me = MException('NCTOOLBOX:ncvariable:somedata', ...
                            ['The data size of the coordinate variable,' ...
                            name ', does not fit the size of ' obj.name]);
                        me.throw;
                    end
                    
                else
                    %% case: variable is coordiantes. Look for size
                    % TODO this is a lame implementation.
                    dim = find(s == vs(1), 1);
                    if ~isempty(dim)
                        for j = 2:length(vs)
                            if vs(j) ~= s(dim + j - 1)
                                me = MException('NCTOOLBOX:ncvariable:somedata', ...
                                    ['The data size of the coordinate variable,' ...
                                    name ', does not fit the size of ' obj.name]);
                                me.throw;
                            end
                        end
                        k = dim:dim + length(vs) - 1;
                        vFirst = first(k);
                        vLast = last(k);
                        vStride = stride(k);
                    end
                end
                
                n = max(size(obj.size(v.getName)));
                ranges = java.util.ArrayList(n);
                for j = 1:n
                    ranges.add(ucar.ma2.Range(vFirst(j) - 1, vLast(j) - 1, vStride(j)));
                end
                
                array = v.read(ranges);
                values = array.copyToNDJavaArray();
                gr.(coordinate_vars{i}) = values;
                clear values
                
            end
            %         else
            % not needed i dont think...
            
            
            
            % Get connectivity
            v = obj.netcdfugrid.getNetcdfFile().findVariable(connectivity_var);
            array = v.read();
            values = array.copyToNDJavaArray();
            gr.connectivity = values;
            
        end % end grid
        
        function vs = size(obj, varName)
            % NCUGRID.SIZE - Return the dimensional sizes of the variable of interest.
            %
            % Usage:
            %              d = nc.size(variableName)
            %
            v = obj.dataset.netcdf.findVariable(varName);
            vs = v.getShape;
            vs = vs';
        end
        
        function d = timewindowij(src, varargin)
            % NCGEOVARIABLE.TIMEWINDOWIJ - Function to get indices from start and stop times for sub-
            % setting. TODO: There must be a better/fast way to do this using the java library.
            % Useage: >> timestruct = nc.timewindowij([2004 1 1 0 0 0], [2005 12 31 0 0 0]);
            %              >> timestruct = nc.timewindowij(731947, 732677);
            % Result: time.time = time data
            %            time.index = time indices
            s = src.size(varargin{1});
            first = ones(1, length(s));
            last = s;
            stride = first;
            g = src.grid_interop(first, last, stride);
            
            if isfield(g, 'time') % are any of the fields recognized as time explictly
                starttime = datenum(varargin{2});
                stoptime = datenum(varargin{3});
                if isempty(starttime)
                    starttime = g.time(1);
                end
                if isempty(stoptime)
                    stoptime = g.time(end);
                end
                
                t_index1 = g.time >= starttime;
                t_index2 = g.time <= stoptime;
                d.index = find(t_index1==t_index2);
                d.time = g.time(d.index);
            else
                me = MException(['NCTOOLBOX:ncgeovariable:timewindowij'], ...
                    'No grid variable returned as time.');
                me.throw;
            end
        end % end timewindowij
        
        
    end % end methods
    
end % end class