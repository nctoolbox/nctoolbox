% POINT Extention of dataset object class for point based datasets -
% Constructor: point(indextype, lat, lon, z_index); which specifies a
% lat/lon or lat/lon index(which in the case of lat/lon index, should be
% the same number since at this point this class assumes that lat and lon
% are the same length)
%
% Alexander Crosby, Applied Science Associates 2010
% 
classdef ncpoint < handle
  % should it be able to handle just points or also collections/grids...?  
  properties (SetAccess = private)
    dataset
%     nativeclass
    
    
  end
  
  methods
    function obj = ncpoint(nc)
%       %Matlab njtbx point class constructor
%       if nargin < 5
%         if nargin < 4
%           indextype = 'index';
%           x = 1;
%           y = 1;
%         end
%         z = ':';
%       end
%       
%       obj.type = indextype;
%       
%       switch obj.type
%         case 'coords'
%           obj.where.lon = x; % I am trying to loosely preserve the 
%                              % x --> lon convention
%           obj.where.lat = y; % y --> lat convention
%           obj.where.z = z; % z is only index right now
%           
%         case 'index'
%           obj.where.lonind = x;
%           obj.where.latind = y;
%           obj.where.z = z; % z is only index right now
%       end
       
      if ischar(nc)
        obj.dataset = ncdataset(src);  % nc is a string URL/File
      elseif isa(nc, 'ncdataset')
        obj.dataset = nc;             % nc is an ncdataset
%       elseif isa(nc, 'ncvariable')
%         obj.dataset = nc;            % nc is an ncvariable
        
      else
        ex = MException('POINT:ncobject', 'Incompatable usage, no reference to local netcdf file, opendap url, ncdataset or ncvariable');
        ex.throw;
      end
        
    end
      
    function ts = timeseries(obj, varargin) % need to add where index to the data call here
    % Usage: 
    % ts = ncpoint.timeseries;  % struct ts, ALL vars ALL time
    % ts = ncpoint.timeseries({variable(s)});  %select vars ALL time
    % ts = ncpoint.timeseries(starttime, stoptime);  %ALL vars btwn
    % start/stop
    % ts = ncpoint.timeseries(starttime, stoptime, {variable(s)});
    
      
      
    % Need to spec the format of start and stop time...
%       starttime = varargin{1};
%       stoptime = varargin{2};
%       variables = varargin{3};
      if nargin < 4 % nargin 0, 1 or 2
        if nargin > 1 % nargin 1 or 2
          if nargin < 3 % nargin of 1
            variables = varargin{1};
            % Take ALL TIME from 1:end
            ts.time = obj.dataset.time('time'); % this is a bad assumption
            % take the LIST of variables
            if iscell(variables)
              for i = 1:length(variables)
                ts.(variables{i}) = obj.dataset.data(variables{i});
              end
            else
              ts.(variables) = obj.dataset.data(variables);
            end
            
          else % nargin of 2
            starttime = datenum(varargin{1});
            stoptime = datenum(varargin{2});
            % Take time series from starttime:stoptime
            t_converted = obj.dataset.time('time'); % this is a bad assumption
            t_index1 = t_converted > starttime;
            t_index2 = t_converted < stoptime;
            t_index = find(t_index1==t_index2);
            ts.time = t_converted(t_index);
            % AND take ALL variables
            vars = obj.dataset.variables;
            for i=1:length(vars)
              dat = obj.dataset.data(vars{i}); % interim notation
              if length(dat)==length(t_converted)
                ts.(vars{i}) = [];
                ts.(vars{i}) = dat(t_index);
              end
              clear dat
            end
          end
          
        else % nargin of 0
          % Take time series from ALL
          % Convert from [yyyy mm dd hh MM ss] to the format specified by the
          % units attribute of time variable
          ts.time = obj.dataset.time('time'); % this is a bad assumption
          % AND---> ALL variables
          vars = obj.dataset.variables;
          for i=1:length(vars)
            dat = obj.dataset.data(vars{i}); % interim notation
            if length(dat)==length(ts.time)
              ts.(vars{i}) = [];
              ts.(vars{i}) = dat;
            end
            clear dat
          end
        end
      else % nargin of 3
        % Take time series froms starttime:stoptime
        starttime = datenum(varargin{1});
        stoptime = datenum(varargin{2});
        variables = varargin{3};
        
        t_converted = obj.dataset.time('time'); % this is a bad assumption
        t_index1 = t_converted > starttime;
        t_index2 = t_converted < stoptime;
        t_index = find(t_index1==t_index2);
        ts.time = t_converted(t_index);
        % Also take the LIST of variables
        if iscell(variables)
          for i = 1:length(variables)
            dat = obj.dataset.data(variables{i});
            ts.(variables{i}) = dat(t_index);
            clear dat
          end
        else
          dat = obj.dataset.data(variables);
          ts.(variables) = dat(t_index);
          clear dat
        end
      end
      % configure output of time and variable values
      % use simple matrix construct: ts(:,1)=time and ts(:,2:end)=values
    end
    
    function s = size(obj, variable, lat, lon)
    % Find size of time and z at any point for reference variable.
    % Usage:
    % s = ncpoint.size(variable)
      
      % find the grid/time variable, find the length at (i,j)
      dat = obj.dataset.data(variable);
      [~ , ~, ind] = obj.nearto(lat, lon, variable);
      
      if length(size(dat)) > 3
        [s1 s2 s3] = size(dat(:,:,ind));
      else 
        [s1 s2 s3] = size(dat(:,ind));
      end
      s = [s1 s2];
      clear dat
    end
    
    function b = bbox(obj, east_min, north_min, east_max, north_max, variable)
    % Get indicies of stations/points that are inside of the box specified.
    % The assumption is that station, lat, lon are all the same length.
    % Usage:
    % ind = ncpoint.bbox(-95, 28, -90, 32)
      
      % Separate from the grid based bbox in ncvariable.
      % Use dataset.grid, find index of lat/lon pairs inside of box use the
      % index to grab the stations that are in the bbox with data. Will
      % have to pull down the entire dataset before I do this in matlab,
      % maybe better way in Java? ----> or not if i index with
      % ncvariable.dataset.data without subsref (hopefully...)
      v = obj.dataset.variable(variable);
      a = v.axes;
      g = v.grid; % need to figure out a way to search
                                      % for variables that arn't of the
                                      % grid
      lat_name = char(a(end-1));
      lon_name = char(a(end));
      indlat = (g.(lat_name) <= north_max)&(g.(lat_name) >= north_min);
      indlon = (g.(lon_name) <= east_max)&(g.(lon_name) >= east_min);
      
      b = find(indlat&indlon);

%       indlat_start = min(indlat);
% 
%       indlat_end = max(indlat);

%       b.data = obj.data(...
%         [1 indlat_start indlon_start],...
%         [nums(1) indlat_end,indlon_end],...
%         [1 1 1]...
%         );
%       % if nargout > 1
%       b.grid = obj.grid(...
%         [1 indlat_start indlon_start],...
%         [nums(1) indlat_end,indlon_end],...
%         [1 1 1]...
%         );
      
    end % not compatible with grids, assumes lat/lon are the same length
    
    function nt = nearto(obj, nlat, nlon, variable) % assumes lat/lon are the same length
    % Method to set/reset the ncpoint.where property of the ncpoint object
    % to find the index and corresponding lat/lon values for the point
    % closest to the the location described by nlon and nlat.
    % Usage:
    % ~ = ncpoint.nearto(42, -70, reference_variable);
    % [lat lon ind] = ncpoint.nearto(42, -70, var);
      
      % get grid (ie lats lons and calculate the distance to each pair
      g = obj.dataset.grid(variable);
      if isfield(g, 'lat')
          lat = g.lat;
          lon = g.lon;
        
      elseif isfield(g, 'Y')
          lat = g.Y;
          lon = g.X;
        
      elseif isfield(g, 'y')
          lat = g.y;
          lon = g.x;
       
      else
        error('bbox:dimensions','No dimension handle found, specifically latitude')
      end
      
      xdist = abs(lon-nlon);
      ydist = abs(lat-nlat);
      hy_dist = sqrt((xdist.^2)+(ydist.^2)); % should be a vector of distances to each station
      closest = find(hy_dist==min(hy_dist)); 
      lon = lon(closest);
      lat = lat(closest);
      
      nt = [lat lon closest];
    end
    
  end
  
  
end