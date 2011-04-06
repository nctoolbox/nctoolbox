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
        
        function ig = interopgrid(src, first, last, stride) % layer subsref for matlab indexing
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
          
        end % interopgrid end
        
        function tw = timewindow(src, starttime, stoptime)
          d = src.timewindowij(starttime, stoptime);
          tw = d.time;
        end
        
        %% These functions would rather output multiple outputs instead of struct, must reconcile
        %     with the subsref in either ncgeovariable or ncvariable. Wait, why, then, does geoij work???
        function d = timewindowij(src, starttime, stoptime)
          s = src.size;
          first = ones(1, length(s));
          last = s;
          stride = first;
          g = src.interopgrid(first, last, stride);
          
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
        
        function d = timegeosubset(src, starttime, stoptime, zmin, zmax, east_min, north_min, ...
            east_max, north_max, stride)
          nums = src.size;
            
            [indstart_r indend_r indstart_c indend_c] = src.geoij(east_min, north_min, east_max, north_max);
            t = src.timewindowij(starttime, stoptime);
            
            if length(nums) < 3
              me = MException(['NCTOOLBOX:' mfilename ':geosubset'], ...
                ['Expected data of ', obj.name, ' to be at least rank 3.']);
              me.throw;
            elseif length(nums) < 4
              first = [min(t.index) indstart_r indstart_c];
              last = [max(t.index) indend_r indend_c];
            elseif length(nums) < 5
              first = [min(t.index) zmin indstart_r indstart_c];
              last = [max(t.index) zmax indend_r indend_c];
            else
              me = MException(['NCTOOLBOX:' mfilename ':geosubset'], ...
                ['Expected data of ', obj.name, ' to be less than rank 5.']);
              me.throw;
            end
            d.data = src.data(first, last, stride);
            d.grid = src.interopgrid(first, last, stride);
            
        end % end of timegeosubset
        
        function d = geosubset(obj, tmin_i, tmax_i, zmin_i, zmax_i, east_min,...
                north_min, east_max, north_max, stride)
            % GEOVARIABLE.GEOSUBSET
            %
            % For use with nj_tbx/nctoolbox to return data based on geographic extents.
            % Use:
            % data = variable.geosubset(1, 1000, 2, 5, -71.5, 39.5, -65, 46, [1 1 1 1])
            %                         %[mintime_ind, maxtime_ind, minZ_ind, maxZ_ind, mineast, minnorth, maxeast, maxnorth, [stride]]%
            %
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
            
            [indstart_r indend_r indstart_c indend_c] = obj.geoij(east_min, north_min, east_max, north_max);
            
            if length(nums) < 2
              me = MException(['NCTOOLBOX:' mfilename ':geosubset'], ...
                ['Expected data of ', obj.name, ' to be at least rank 2.']);
              me.throw;
            elseif length(nums) < 3
              first = [indstart_r indstart_c];
              last = [indend_r indend_c];
            elseif length(nums) < 4
              first = [tmin_i indstart_r indstart_c];
              last = [tmax_i indend_r indend_c];
            elseif length(nums) < 5
              first = [tmin_i zmin_i indstart_r indstart_c];
              last = [tmax_i zmax_i indend_r indend_c];
            else
              me = MException(['NCTOOLBOX:' mfilename ':geosubset'], ...
                ['Expected data of ', obj.name, ' to be less than rank 5.']);
              me.throw;
              
            end
            d.data = obj.data(first, last, stride);
            d.grid = obj.interopgrid(first, last, stride);
            
        end % end of geosubset
       %% 
        function [indstart_r indend_r indstart_c indend_c] =...
                geoij(obj, east_min, north_min, east_max, north_max)
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
            g = obj.interopgrid(first, last, stride);
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
            
            indstart_c = min(ind_c);
            indend_c = max(ind_c);
            indstart_r = min(ind_r);
            indend_r = max(ind_r);

        end
    end % methods end
    
%     methods (Access = protected)
%      
%     end % protected methods end
end % class end