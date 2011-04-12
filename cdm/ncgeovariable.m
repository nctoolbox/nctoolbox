% NCGEOVARIABLE Provide advanced access to variables and their related
% dimensions.
%
% NCGEOVARIABLE is used to retrieve data for a given variable as well as the
% variables associated coordinate dimensions.
% 
%
% Example of use:
%  ds = cfdataset('http://dods.mbari.org/cgi-bin/nph-nc/data/ssdsdata/deployments/m1/200810/OS_M1_20081008_TS.nc');
%  v = ds.variable('TEMP');
%  t = v.data([1 1 1 1], [100 5 1 1]);
%  % Look at properties
%  v.name
%  v.axes

classdef ncgeovariable < ncvariable
    
    properties (SetAccess = private)
%         dataset          % ncdataset instance
    end
    
    properties (Dependent = true)
%         name            % The string variable name that this object represents
%         attributes
    end
    
    properties (SetAccess = private, GetAccess = protected)
%         variable        % ucar.nc2.Variable instance. Represents the data
%         axesVariables    % ucar.nc2.Variable instance. Represents the data.

    end
    
    methods
        
        %%
        function obj = ncgeovariable(src, variableName, axesVariableNames)
            % NCGEOVARIABLE.NCGEOVARIABLE  Constructor.
            %
            % Use as:
            %    v = ncvariable(src, variableName)
            %    v = ncvariable(src, variableName, axesVariableNames)
            %
            obj = obj@ncvariable(src, variableName, axesVariableNames);
            
            
        end % ncgeovariable end
        
        function ig = grid_interop(src, first, last, stride) % layer subsref for matlab indexing
          g = src.grid(first, last, stride);
          names = fieldnames(g);
          
          for i = 1:length(names); % loop through fields returned by grid
            tempname = names{i};
            javaaxisvar  =   src.dataset.netcdf.findVariable(tempname);
            type = char(javaaxisvar.getAxisType());
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
          end % end loop through field names
          
        end % grid_interop end
        
        function tw = timewindow(src, starttime, stoptime)
          d = src.timewindowij(starttime, stoptime);
          tw = d.time;
        end

        
        %% These functions would rather output multiple outputs instead of struct, must reconcile
        %     with the subsref in either ncgeovariable or ncvariable. Wait, why, then, does geoij work???
        function d = timewindowij(src, starttime, stoptime)
          % NCGEOVARIABLE.TIMEWINDOWIJ - Function to get indices from start and stop times for sub-
          % setting. TODO: There must be a better/fast way to do this using the java library.
          s = src.size;
          first = ones(1, length(s));
          last = s;
          stride = first;
          g = src.grid_interop(first, last, stride);
          
          if isfield(g, 'time') % are any of the fields recognized as time explictly
            starttime = datenum(starttime);
            stoptime = datenum(stoptime);
            if isempty(starttime)
              starttime = g.time(1);
            end
            if isempty(stoptime)
              stoptime = g.time(end);
            end
            
            t_index1 = g.time > starttime;
            t_index2 = g.time < stoptime;
            d.index = find(t_index1==t_index2);
            d.time = g.time(d.index);
          else
            me = MException(['NCTOOLBOX:' mfilename ':timewindowij'], ...
                'No grid variable returned as time.');
            me.throw;
          end
        end % end timewindowij
        
        function d = timegeosubset(src, struct)
          nums = src.size;
            
            [indstart_r indend_r indstart_c indend_c] = src.geoij(struct);
            t = src.timewindowij(struct.time{1}, struct.time{2});
            
            if length(nums) < 3
              me = MException(['NCTOOLBOX:' mfilename ':geosubset'], ...
                ['Expected data of ', obj.name, ' to be at least rank 3.']);
              me.throw;
            elseif length(nums) < 4
              ax = obj.grid([1 1 1],[1 1 1],[1 1 1]);
              if isfield(ax, 'time')
                first = [min(t.index) indstart_r indstart_c];
                last = [max(t.index) indend_r indend_c];
                stride = [struct.t_stride struct.xy_stride(2) struct.xy_stride(1)];
              else
                me = MException(['NCTOOLBOX:' mfilename ':geosubset'], ...
                  'Expected either a coordinate variable acknowleged as time.');
                me.throw;
              end
            elseif length(nums) < 5
              first = [min(t.index) struct.z_index{1} indstart_r indstart_c];
              last = [max(t.index) struct.z_index{2} indend_r indend_c];
              stride = [struct.t_stride struct.z_stride struct.xy_stride(2) struct.xy_stride(1)];
            else
              me = MException(['NCTOOLBOX:' mfilename ':geosubset'], ...
                ['Expected data of ', obj.name, ' to be less than rank 5.']);
              me.throw;
            end
            d.data = src.data(first, last, stride);
            d.grid = src.grid_interop(first, last, stride);
            
        end % end of timegeosubset
        
        function d = geosubset(obj, struct)
            % GEOVARIABLE.GEOSUBSET
            %
            % For use with nj_tbx/nctoolbox to return data based on geographic extents. 
            % Use: data = variable.geosubset(struct); %Where struct is the kind of structure produced by 
            %         geosubset_struct.m.
            %                         
            %
            %
            % TODO: add stride arguments and catches for points and stations
            % because this logic won't work with them.
            % Alexander Crosby, Applied Science Associates
            %
            nums = obj.size;
            
            [indstart_r indend_r indstart_c indend_c] = obj.geoij(struct);
            
            if numel(struct.time{1}) > 1 % check to see if someone used str or datevec by accident
              me = MException(['NCTOOLBOX:' mfilename ':geosubset'], ...
                'Expected min time to be an index/integer.');
              me.throw;
            else
              tmin_i = struct.time{1};
            end
            
            if numel(struct.time{2}) > 1 % check to see if someone used str or datevec by accident
              me = MException(['NCTOOLBOX:' mfilename ':geosubset'], ...
                'Expected max time to be an index/integer.');
              me.throw;
            else
              tmax_i = struct.time{2};
            end
            
            if length(nums) < 2
              me = MException(['NCTOOLBOX:' mfilename ':geosubset'], ...
                ['Expected data of ', obj.name, ' to be at least rank 2.']);
              me.throw;
            elseif length(nums) < 3
              first = [indstart_r indstart_c];
              last = [indend_r indend_c];
              stride = [struct.xy_stride(2) struct.xy_stride(1)];
            elseif length(nums) < 4
              ax = obj.grid([1 1 1],[1 1 1],[1 1 1]);
              if isfield(ax, 'time')
                first = [tmin_i indstart_r indstart_c];
                last = [tmax_i indend_r indend_c];
                stride = [struct.t_stride struct.xy_stride(2) struct.xy_stride(1)];
              elseif isfield(ax, 'z')
                first = [struct.z_index{1} indstart_r indstart_c];
                last = [struct.z_index{2} indend_r indend_c];
                stride = [struct.z_stride struct.xy_stride(2) struct.xy_stride(1)];
              else
                me = MException(['NCTOOLBOX:' mfilename ':geosubset'], ...
                  'Expected either a coordinate variable acknowleged as time or as z.');
                me.throw;
              end
            elseif length(nums) < 5
              first = [tmin_i struct.z_index{1} indstart_r indstart_c];
              last = [tmax_i struct.z_index{2} indend_r indend_c];
              stride = [struct.t_stride struct.z_stride struct.xy_stride(2) struct.xy_stride(1)];
            else
              me = MException(['NCTOOLBOX:' mfilename ':geosubset'], ...
                ['Expected data of ', obj.name, ' to be less than rank 5.']);
              me.throw;
              
            end
            
            % Get the corresponding data and interop grid...
            d.data = obj.data(first, last, stride);
            d.grid = obj.grid_interop(first, last, stride);
            
        end % end of geosubset
       %% 
        function [indstart_r indend_r indstart_c indend_c] =...
                geoij(obj, struct)
            % GEOVARIABLE.GEOIJ
            %
            % For use with nj_tbx/nctoolbox to return data based on geographic extents.
            %
            % This code relys on coards conventions of coodinate order using:
            % [time, z, lat, lon]
            %
            % TODO: add stride arguments and catches for points and stations
            % because this logic won't work with them.
            % Alexander Crosby, Applied Science Associates
            %
            s = obj.size;
            first = ones(1, length(s));
            last = s;
            stride = first;
            g = obj.grid_interop(first, last, stride);
            %           h = 0;
            
            %Unpack geosubset_structure
            north_max = struct.lat(2);
            north_min = struct.lat(1);
            east_max = struct.lon(2);
            east_min = struct.lon(1);
            
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
            
            indstart_c = min(ind_c);
            indend_c = max(ind_c);
            indstart_r = min(ind_r);
            indend_r = max(ind_r);

        end
        
        function sref = subsref(obj,s)
            switch s(1).type
                % Use the built-in subsref for dot notation
                case '.'
                    switch s(1).subs
                      case 'grid_interop'
                      switch length(s)
                        case 1
                          sref = obj;
                        case 2
                          nums = obj.size;
                          [first last stride] = parseIndices(s(2).subs, double(nums));
                          sref = obj.grid_interop(first, last, stride);
                      end
                      case 'data'
                            nums = size(obj);
                            if ~isempty(nums)
                                switch length(s)
                                    case 1
                                        sref = obj;
                                    case 2
                                        [first last stride] = parseIndices(s(2).subs, double(nums));
                                        sref = obj.data(first, last, stride);
                                end
                                
                            else
                                sref = obj.data;
                                warning(['NCTOOLBOX:' mfilename ':subsref'], ...
                                    ['Variable "' name '" has no netcdf dimension associated with it. Errors may result from non CF compliant files.'])
                            end
                        case 'grid'
                            nums = size(obj);
                            if ~isempty(nums)
                                switch length(s)
                                    case 1
                                        sref = obj;
                                    case 2
                                        [first last stride] = parseIndices(s(2).subs, double(nums));
                                        sref = obj.grid(first, last, stride);
                                end

                            else
                                warning(['NCTOOLBOX:' mfilename ':subsref'], ...
                                    ['Variable "' name '" has no netcdf dimension associated with it. Errors may result from non CF compliant files.'])
                            end
                      otherwise
                        sref = builtin('subsref',obj,s);
                    end
                case '()'
                    warning(['NCTOOLBOX:' mfilename ':subsref'], ...
                        'Not a supported subscripted reference, "()" are not permitted to call variable object methods');
                case '{}'
                    warning(['NCTOOLBOX:' mfilename ':subsref'], ...
                        'Not a supported subscripted reference, "{}" are not permitted to call variable object methods');
            end
        end
    end % methods end
    
%     methods (Access = protected)
%      
%     end % protected methods end
end % class end