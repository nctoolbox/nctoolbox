% NCVARIABLE Provide advanced access to variables and their related
% dimensions.
%
% NCVARIABLE is used to retrieve data for a given variable as well as the 
% variables associated coordinate dimensions. Normally, you would retrive
% it using CFDATASET.VARIABLE 
%
% Example of use:
%  ds = cfdataset('http://dods.mbari.org/cgi-bin/nph-nc/data/ssdsdata/deployments/m1/200810/OS_M1_20081008_TS.nc');
%  v = ds.variable('TEMP');
%  t = v.data([1 1 1 1], [100 5 1 1]);
%  % Look at properties
%  v.name
%  v.axes
%
%
% See also CFDATASET, SIZE, DATA
classdef ncvariable < handle
    
    properties (SetAccess = private)
        dataset          % ncdataset instance
    end
    
    properties (Dependent = true)
        name            % The string variable name that this object represents
        axes
    end
    
    properties (SetAccess = private, GetAccess = protected)
        variable        % ucar.nc2.Variable instance. Represents the data
        axesVariables    % ucar.nc2.Variable instance. Represents the data.
        attributes %added acrosby
    end
    
    methods
        
        %%
        function obj = ncvariable(src, variableName, axesVariableNames)
            % NCVARIABLE.NCVARIABLE  Constructor.
            %
            % Use as:
            %    v = ncvariable(src, variableName)
            %    v = ncvariable(src, variableName, axesVariableNames)
            %
            if ischar(src)
                obj.dataset = ncdataset(src);  % src is a string URL/File
            elseif isa(src, 'ncdataset')
                obj.dataset = src;             % src is a ncdataset
            else
                ex = MException('NCVARIABLE:ncvariable', 'Invalid dataset was specified');
                ex.throw;
            end
            
            obj.variable = obj.dataset.netcdf.findVariable(variableName);
            obj.attributes = obj.dataset.attributes(obj.name); %acrosby
            
            if nargin == 3
                obj.axesVariables = cell(size(axesVariableNames));
                for i = 1:length(axesVariableNames)
                    obj.axesVariables{i} = obj.dataset.netcdf.findVariable(axesVariableNames{i});
                end
            else
                obj.axesVariables = {};
            end
            
            
        end
        
        %%
        function a = get.axes(obj)
            % NCVARIABLE.AXES Returns the names of the coordinate axes 
            a = cell(size(obj.axesVariables));
            for i = 1:length(obj.axesVariables)
                name = char(obj.axesVariables{i}.getName());
                a{i} = name;
            end
        end
        
        %%
        function v = get.name(obj)
            % NCVARIABLE.NAME Provides dynamic access to the underlying
            % netcdf datasets variable name
            v = char(obj.variable.getName());
        end
        
        %%
        function n = size(obj)
            % NCVARIABLE.SIZE returns the size of the variable, including
            % it's singleton dimensions
            n = obj.dataset.size(obj.name);
        end
        
        %%
        function d = data(obj, first, last, stride)
            % NCVARIABLE.DATA Retrieve all or a subset of the data for the
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
            %
            % Returns:
            %   The data is returned as a structure containing the actual data for the variable 
            %   of interest as well as each coordinate variable
            %
            % Example:
            %
            %   ds = cfdataset('http://dods.mbari.org/cgi-bin/nph-nc/data/ssdsdata/deployments/m1/200810/OS_M1_20081008_TS.nc');
            %   v = ds.variable('TEMP');
            %   t = v.data([1 1 1 1], [10 2 1 1]);
            %

            if (nargin == 1)
              if isempty(obj.size) % check if variable obj is dimensionless (not cf compliant.....)
                try
                  d = alldata(obj,1);
                catch me
                  me.throw();
                end
              else
                d = alldata(obj, 1);
              end
            else
                s = obj.size;
                
                % Fill in missing arguments
                % default stride is 1
                if (nargin < 4)
                    stride = ones(1, length(s));
                end
                
                % Default last is the end
                if (nargin < 3)
                    last = s;
                end
                
                d = somedata(obj, 1, first, last, stride);
            end
        end
        
        function g = grid(obj, first, last, stride)
            % NCVARIABLE.GRID Retrieve all or a subset of the coordinate 
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
            %   The data is returned as a structure containing the actual data for the variable 
            %   of interest as well as each coordinate variable
            %
            % Example:
            %
            %   ds = cfdataset('http://dods.mbari.org/cgi-bin/nph-nc/data/ssdsdata/deployments/m1/200810/OS_M1_20081008_TS.nc');
            %   v = ds.variable('TEMP');
            %   t = v.data([1 1 1 1], [10 2 1 1]);
            %
            
            if (nargin == 1)
                g = alldata(obj, 0);
            else % all this stuff needs logic work to function as expected
                s = obj.size;
                
                % Fill in missing arguments
                % default stride is 1
                if (nargin < 4)
                    stride = ones(1, length(s));
                end
                
                % Default last is the end
                if (nargin < 3)
                    last = s;
                end
                
                g = somedata(obj, 0, first, last, stride);
            end
        end
        
        function t = gettime(obj, varargin)
          g = obj.grid;
          if isfield(g, 'time')
            t_converted = g.time; % this is a bad assumption
            if nargin > 2
              t_index1 = t_converted > varargin{2};
              t_index2 = t_converted < varargin{3};
              t_index = find(t_index1==t_index2);
              t = t_converted(t_index);
            else
              t = t_converted;
            end
          else
            error('gettime:ncvariable','No variable that netcdf-java recognizes as time.');
          end
        end
        
        %%
        function bb = bbox(obj, tmin_i, tmax_i, zmin_i, zmax_i, east_min,...
            north_min, east_max, north_max, stride)
          % BBOX
          %
          % For use with nj_tbx/nctoolbox to return data based on geographic extents.
          % Use:
          % data = variable.bbox(-71.5, 39.5, -65, 46)
          %                         %[east north east north]%
          %
          % This code relys on coards conventions of coodinate order using:
          % [time, z, lat, lon]
          %
          % TODO: add stride arguments and catches for points and stations
          % because this logic won't work with them.
          % Alexander Crosby, Applied Science Associates
          %
%           g = obj.grid;
%           h = 0;
%           a = obj.axes;
%           [lat_name] = char(a(end-1));
%           [lon_name] = char(a(end));
          nums = obj.size;
          
          [indstart_r indend_r indstart_c indend_c] = obj.bboxij(east_min, north_min, east_max, north_max);
          
          if length(nums) < 4

            tstart = first(1);
            tend = last(1);
            zstart = first(2);
            zend = last(2);
            first = [1 indstart_r indstart_c];
            last = [1 indend_r indend_c];
%             stride = [1 1 1];
            
          else
            first = [tmin_i zmin_i indstart_r indstart_c];
            last = [tmax_i zmax_i indend_r indend_c];
%             stride = [1 1 1 1];
            
          end
          bb.data = obj.data(first, last, stride);
          bb.grid = obj.grid(first, last, stride);
          
        end 
        
        function [indstart_r indend_r indstart_c indend_c] =...
            bboxij(obj, east_min, north_min, east_max, north_max)
          % BBOXIJ
          %
          % For use with nj_tbx/nctoolbox to return data based on geographic extents.
          % Use:
          % data = variable.bbox(-71.5 39.5 -65 46)
          %                         %[east north east north]%
          %
          % This code relys on coards conventions of coodinate order using:
          % [time, z, lat, lon]
          %
          % TODO: add stride arguments and catches for points and stations
          % because this logic won't work with them.
          % Alexander Crosby, Applied Science Associates
          %
          g = obj.grid;
%           h = 0;
          
          

            if ~isvector(g.lat)
              [indlat_l1] = ((g.lat <= north_max)); %2d
              [indlat_l2] = ((g.lat >= north_min)); %2d
              [indlat_r indlat_c] = find((indlat_l1&indlat_l2)); % 1d each
              indlon_l1 = zeros(size(indlat_l1));
              indlon_l2 = zeros(size(indlat_l1));
              for i = 1:length(indlat_r)
                if g.lon(indlat_r(i), indlat_c(i)) <= east_max;
                  indlon_l1(indlat_r(i), indlat_c(i)) = true;                
                end
                if g.lon(indlat_r(i), indlat_c(i)) >= east_min;
                  indlon_l2(indlat_r(i), indlat_c(i)) = true;                
                end
              end
              [ind_r, ind_c] = find((indlon_l1&indlon_l2));
              h=1;
            else
              indlat1 = (g.lat <= north_max);
              indlat2 = (g.lat >= north_min);
              indlat = find(indlat1&indlat2);
              indlon1 = (g.lon <= east_max);
              indlon2 = (g.lon >= east_min);
              indlon = find(indlon1&indlon2);
              
            end
            
%             if h==1 %if lat, lon are 2-D  %% doesn't matter?
            indstart_c = min(ind_c);
            indend_c = max(ind_c);
            indstart_r = min(ind_r);
            indend_r = max(ind_r);
%             nums = obj.size;
%             if length(nums) < 4
%               first = [1 indstart_r indstart_c];
%               last = [1 indend_r indend_c];
%               stride = [1 1 1];
% 
%             else
%               first = [1 1 indstart_r indstart_c];
%               last = [1 1 indend_r indend_c];
%               stride = [1 1 1 1];
% 
%             end
%             bb.data = obj.data(first, last, stride);
%             bb.grid = obj.grid(first, last, stride);
            
%           else %if lat, lon are 1-D
%             indlat_start = min(indlat);
%             indlon_start = min(indlon);
%             indlat_end = max(indlat);
%             indlon_end = max(indlon);
% %             nums = size(obj);
% 
% %             if length(nums) < 4
% %             bb.data = obj.data(...
% %               [1 indlat_start indlon_start],...
% %               [nums(1) indlat_end,indlon_end],...
% %               [1 1 1]...
% %               );
% %             %         if nargout > 1
% %             bb.grid = obj.grid(...
% %               [1 indlat_start indlon_start],...
% %               [nums(1) indlat_end,indlon_end],...
% %               [1 1 1]...
% %               );
% %             else
% %               bb.data = obj.data(...
% %               [1 1 indlat_start indlon_start],...
% %               [nums(1) nums(2) indlat_end,indlon_end],...
% %               [1 1 1 1]...
% %               );
% %             %         if nargout > 1
% %               bb.grid = obj.grid(...
% %               [1 1 indlat_start indlon_start],...
% %               [nums(1) nums(2) indlat_end,indlon_end],...
% %               [1 1 1 1]...
% %               );
% %             end
          end
        

function sref = subsref(obj,s)
  switch s(1).type
    % Use the built-in subsref for dot notation
    case '.'
      switch s(1).subs
        case 'data'
          nums = size(obj);
          if ~isempty(nums)
            switch length(s)
              case 1
                first = ones(1,length(nums));
                last = double(nums);
                stride = first;
              case 2
%                 switch length(nums)
%                   case 1
                    [first last stride] = parseIndices(s(2).subs, double(nums));

            end
            sref = obj.data(first, last, stride);
          else
            sref = obj.data;
            disp('Warning: Variable has no netcdf dimension associated with it. Errors may result from non CF compliant files.')
          end
        case 'grid'
          nums = size(obj);
          if ~isempty(nums)
            switch length(s)
              case 1
                first = ones(1,length(nums));
                last = double(nums);
                stride = first;
              case 2
%                 switch length(nums)
%                   case 1
                    [first last stride] = parseIndices(s(2).subs, double(nums));
            end
            sref = obj.grid(first, last, stride);
          else
            disp('No dimensions associated with variable');
          end
        otherwise
          sref = builtin('subsref',obj,s);
      end
    case '()'
      if length(s)<2
        % Note that obj.Data is passed to subsref
        sref = builtin('subsref',obj.Data,s);
      else
        sref = builtin('subsref',obj,s);
      end
      % No support for indexing using '{}'
    case '{}'
      error('variable_object:subsref',...
        'Not a supported subscripted reference, "{}" are not permitted to call variable object methods')
  end
end


    end
    
    methods (Access = protected)
        
       
        %%
        function data = alldata(obj, withData)
            
            % ---- Step 2: Add the data
            if withData == 1
                name = char(obj.variable.getName());
                data = obj.dataset.data(name);
            end
            
             if withData == 0
              for i = 1:length(obj.axesVariables)
                name = char(obj.axesVariables{i}.getName());
                data.(name) = obj.dataset.data(name);
              end
            end
        end
        
        %%
        function data = somedata(obj, withData, first, last, stride)
            
            s = obj.dataset.size(obj.name);
            
            % Fill in missing arguments
            % default stride is 1
            if (nargin < 5)
                stride = ones(1, length(s));
            end
            
            % Default last is the end
            if (nargin < 4)
                last = s;
            end
            
            % ---- Step 2: Add the data for the variable of interest
            if withData==1
                name = char(obj.variable.getName());
                data = obj.dataset.data(name, first, last, stride);
            end

            % ---- Step 3: Add the data for each axes variable
            if withData==0
            for i = 1:length(obj.axesVariables)
                name = char(obj.axesVariables{i}.getName());
                type = char(obj.axesVariables{i}.getAxisType());
                % ---- Step 4: figure out how to subset the data properly
                vs = obj.dataset.size(name);
                if (length(vs) == length(s)) 
                    %% case: sizes are the same
                    if isequal(vs, s)
                        vFirst = first;
                        vLast = last;
                        vStride = stride;
                    else
                        me = MException('NCVARIABLE:somedata', ['The data size of the coordinate variable,' ...
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
                        me = MException('NCVARIABLE:somedata', ['The data size of the coordinate variable,' ...
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
                                me = MException('NCVARIABLE:somedata', ['The data size of the coordinate variable,' ...
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
                
                if isempty(type)
                  data.(name) = obj.dataset.data(name, vFirst, vLast, vStride);
                else
                  switch type
                    case 'Height'
                      pos_z = char(obj.axesVariables{i}.getPositive());
                      if strcmp(pos_z, 'POSITIVE_DOWN')
                        tmp = obj.dataset.data(name, vFirst, vLast, vStride);
                        data.z = tmp.*-1; %adjust for positive direction
                      else
                        data.z = obj.dataset.data(name, vFirst, vLast, vStride);
                      end
                      
                      case 'GeoZ'
                      pos_z = char(obj.axesVariables{i}.getPositive());
                      if strcmp(pos_z, 'POSITIVE_DOWN')
                        tmp = obj.dataset.data(name, vFirst, vLast, vStride);
                        data.z = tmp.*-1; %adjust for positive direction
                      else
                        data.z = obj.dataset.data(name, vFirst, vLast, vStride);
                      end
                      
                    case 'Time'
                      tmp = obj.dataset.data(name, vFirst, vLast, vStride);
                      t_converted = obj.dataset.time(name, tmp);
                      data.time = t_converted;
                      
%                     case 'RunTime'
%                       tmp = obj.dataset.data(name, vFirst, vLast, vStride);
%                       t_converted = obj.dataset.time(name, tmp);
%                       data.(type) = t_converted;
                      
                    case 'Lon'
                      tmp = obj.dataset.data(name, vFirst, vLast, vStride);
                      ind = find(tmp > 180);
                      tmp(ind) = tmp(ind)-360;
                      data.lon = tmp;
                    case 'Lat'
                      data.lat = obj.dataset.data(name, vFirst, vLast, vStride);
                    otherwise
                      data.(type) = obj.dataset.data(name, vFirst, vLast, vStride);
                  end
                end
                
            end
            end
        end
    end
    
end