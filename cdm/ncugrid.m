% UGRID Extention of dataset object class for unstructured grid datasets.
% Alexander Crosby, Applied Science Associates 2010
% 
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
        form = Formatter();
        cancelTask = [];
        obj.netcdfugrid = FeatureDatasetFactoryManager.open(FeatureType.UGRID, nc, cancelTask, form);
      
        elseif isa(nc, 'ncdataset')
        obj.dataset = nc;             % src is a ncdataset
        form = Formatter();
        cancelTask = [];
        obj.netcdfugrid = FeatureDatasetFactoryManager.open(FeatureType.UGRID, nc.link, cancelTask, form);
      else
        ex = MException('NCVARIABLE:ncvariable', 'Invalid dataset was specified');
        ex.throw;
      end
      
      %set properties
      u = obj.netcdfugrid;
      obj.meshes = u.getUGrids();
      meshes = obj.meshes;
      obj.numberofgrids = meshes.size();
      for i = 1:(obj.numberofgrids)
        obj.cells(i, 1) = meshes.get(i - 1).getSize();
        obj.nodes(i, 1) = meshes.get(i - 1).getNodeSize();
        obj.edges(i, 1) = meshes.get(i - 1).getEdgeSize();
        obj.faces(i, 1) = meshes.get(i - 1).getFaceSize();
      end
      
      
    end % end constructor
    
    function ss = unstructuredLatLonSubset(obj, varName, struct) % input subset structure
      % NCUGRID.unstructuredLatLonSubset  - Function to subset an unstructured model grid by lat/lon
      % bounding box using subsetting methods in the ugrid-java.
      import ucar.nc2.dt.ugrid.geom.LatLonPolygon2D;
      import ucar.nc2.dt.ugrid.geom.LatLonRectangle2D;
      import ucar.unidata.geoloc.LatLonPoint;
      import ucar.unidata.geoloc.LatLonPointImpl;
      import ucar.unidata.geoloc.LatLonRect;
      % Need to construct lat/lon rect:
      maxlat  = max(struct.lat);
      minlon = min(struct.lon);
      minlat = min(struct.lat);
      maxlon = max(struct.lon);
      bbox = LatLonRect(LatLonPointImpl(maxlat, minlon), LatLonPointImpl(minlat, maxlon));
      
      % Find the mesh that contains the variable of interest.
      mesh_ident = obj.attribute('mesh', varName);
      mesh_ind = str2double(mesh_ident(5:end));
      
      % Get all the meshes, access the correct mesh for the variable, then call subset method on mesh.
      meshes = obj.meshes;
      meshes.get(mesh_ind-1).buildRTree();
      subsat = meshes.get(mesh_ind-1).subset(bbox); % Subsat is also mesh
      
      
      % Get netcdfj variable class and read all (since we already subset)
      %       inds = subsat.getNodeIndexes;
      nodes = subsat.getUniqueNodes();
      % Kyle made a method to do this automatically mesh.getNodeLatLons() and getNodeIndexes()
      for i = 1:nodes.size();
        inds(i,1) = nodes.get(i-1).getDataIndex();
      end
      s = obj.size(varName);
      first = ones(1, length(s));
      last = s;
      stride = first;
      for j = 1:length(inds)
        first(1,end) = inds(j);
        last(1, end) = inds(j);
        ss.data(:,j) = obj.data(varName, first, last, stride);
      end
      % Related coordinates and connectivity 
      ss.grid = obj.getMeSomeCoordinateData(varName, subsat);
      
     
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
        warning('NCTOOLBOX:ncugrid:attribute', 'No key or variable specified.');
      end
    end % end attribute
    
    function gs = getSomeCoordinateData(obj, varName, mesh)
      mesh_ident = obj.attribute('mesh', varName);
      if nargin < 3
        mesh  = obj.meshes.get(str2double(mesh_ident(5:end)-1));
      end

      coordinates = obj.attribute('coordinates', varName);
      coordinate_vars = regexp(coordinates, ' ', 'split');
      
      location = obj.attribute('location', varName);
      connectivity_var = obj.attribute([location, '_connectivity'], mesh_ident);
      
      % Get coordinates (time, lat, lon,...)
      for i = 1:length(coordinate_vars)
        switch coordinate_vars{i}
          case 'time'
            v = mesh.findVariable(coordinate_vars{i});
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
      
      nodes = mesh.getUniqueNodes(); 
      % Kyle made a method to do this automatically mesh.getNodeLatLons() and getNodeIndexes()
      for i = 1:nodes.size();
        point = nodes.get(i-1).getGeoPoint();
        gs.lat(i,1) = point.getLatitude();
        gs.lon(i,1) = point.getLongitude();
        gs.index(i,1) = nodes.get(i-1).getDataIndex();
      end
      % Get connectivity
        v = obj.netcdfugrid.getNetcdfFile().findVariable(connectivity_var);
        array = v.read();
        values = array.copyToNDJavaArray();
        gs.connectivity = values';
        
    end % end grid
    
    function d = data(obj, variable, first, last, stride)
      
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
        d = double(d);
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
      % NCGEOVARIABLE.GRID_INTEROP - Method to get the coordinate variables and their data as a
      % a structure with standardized field names for lat, lon, time, and z. Other coordiante variables
      % that are not recognized by netcdf-java as the previous four types have field names directly
      % taken from their variable names.
      % Useage: >> gridstruct = geovar.grid_interop(1,:,:,1:2:50);
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
      mesh_ident = obj.attribute('mesh', varName);
      if nargin < 3
        mesh  = obj.meshes.get(str2double(mesh_ident(5:end)-1));
      end
      
      coordinates = obj.attribute('coordinates', varName);
      coordinate_vars = regexp(coordinates, ' ', 'split');
      
      location = obj.attribute('location', varName);
      connectivity_var = obj.attribute([location, '_connectivity'], mesh_ident);
      
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
      v = obj.dataset.netcdf.findVariable(varName);
      vs = v.getShape;
      vs = vs';
    end
    
    
%     function e = end(obj, k, n)
% %            n = obj.temp;
%            e = k-1;
%          end % Added to deal with end indexing functionality,
%                                                      % otherwise the indexing arugment is ignored.
%     function sref = subsref(obj,s)
%             switch s(1).type
%                 % Use the built-in subsref for dot notation
%                 case '.'
%                     switch s(1).subs
%                       case 'grid_interop'
%                       switch length(s)
%                         case 1
%                           sref = obj;
% %                           obj.temp = size(s(2).subs{1});
%                         case 2
%                           nums = obj.size(s(2).subs{1});
%                           [first last stride] = indexing(s(2).subs{2:end}, double(nums));
%                           sref = obj.grid_interop(s(2).subs{1}, first, last, stride);
%                       end
%                       case 'data'
%                             
% %                             if ~isempty(nums)
%                                 switch length(s)
%                                     case 1
%                                         sref = obj;
% %                                         obj.temp = size(s(2).subs{1});
%                                     case 2
%                                         nums = size(s(2).subs{1});
%                                         [first last stride] = indexing(s(2).subs{2:end}, double(nums));
%                                         sref = obj.data(s(2).subs{1}, first, last, stride);
%                                 end
%                                 
% %                             else
% %                                 sref = obj.data;
% %                                 warning(['NCTOOLBOX:ncgeovariable:subsref'], ...
% %                                     ['Variable "' name '" has no netcdf dimension associated with it. Errors may result from non CF compliant files.'])
% %                             end
%                         case 'grid'
%                             
% %                             if ~isempty(nums)
%                                 switch length(s)
%                                     case 1
%                                         sref = obj;
% %                                         obj.temp = size(s(2).subs{1});
%                                     case 2
%                                         nums = size(s(2).subs{1});
%                                         [first last stride] = indexing(s(2).subs{2:end}, double(nums));
%                                         sref = obj.grid(s(2).subs{1}, first, last, stride);
%                                 end
% 
% %                             else
% %                                 warning(['NCTOOLBOX:ncgeovariable:subsref'], ...
% %                                     ['Variable "' name '" has no netcdf dimension associated with it. Errors may result from non CF compliant files.'])
% %                             end
%                       otherwise
%                         sref = builtin('subsref',obj,s);
%                     end
%                 case '()'
%                     warning(['NCTOOLBOX:ncgeovariable:subsref'], ...
%                         'Not a supported subscripted reference, "()" are not permitted to call variable object methods');
%                 case '{}'
%                     warning(['NCTOOLBOX:ncgeovariable:subsref'], ...
%                         'Not a supported subscripted reference, "{}" are not permitted to call variable object methods');
%             end
%         end
    
  end % end methods

end % end class