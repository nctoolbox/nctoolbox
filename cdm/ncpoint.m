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
      % NCPOINT.ncpoint - Constructor method for ncpoint class.
      % Usage:
      %            nc = ncpoint('http://testbedapps.sura.org/thredds/dodsC/estuarine_hypoxia/obs/points.nc');
      %            nc = ncpoint(ncdataset);
      if ischar(nc)
        obj.dataset = ncdataset(nc);  % nc is a string URL/File
      elseif isa(nc, 'ncdataset')
        obj.dataset = nc;             % nc is an ncdataset
%       elseif isa(nc, 'ncvariable')
%         obj.dataset = nc;            % nc is an ncvariable
        
      else
        ex = MException('POINT:ncobject', 'Incompatable usage, no reference to local netcdf file, opendap url, ncdataset or ncvariable');
        ex.throw;
      end
        
    end
      
    function ts = timeseries(obj, varargin)
      % NCPOINT.timeseries - Function to subset the data based on a timewindow for all variables based
      % on a given time variable, or a subset of input variables.
      % Usage:
      % ts = ncpoint.timeseries('time');  % struct ts, ALL vars ALL time using time variable 'time'
      % ts = ncpoint.timeseries('time', {variable(s)});  %select vars ALL time using time variable 'time'
      % ts = ncpoint.timeseries('time', starttime, stoptime);  %ALL vars btwn start/stop times with var-
      % iable 'time'
      % ts = ncpoint.timeseries('time', starttime, stoptime, {variable(s)});
      % Need to spec the format of start and stop time...
      tvar = varargin{1};
      if nargin < 5 % nargin 0, 1 or 2
        if nargin > 2 % nargin 1 or 2
          if nargin < 4 % nargin of 1
            variables = varargin{2};
            % Take ALL TIME from 1:end
            
            ts.time = obj.dataset.time(tvar); % this is a bad assumption
            % take the LIST of variables
            if iscell(variables)
              for i = 1:length(variables)
                ts.(variables{i}) = obj.dataset.data(variables{i});
              end
            else
              ts.(variables) = obj.dataset.data(variables);
            end
            
          else % nargin of 2
            starttime = datenum(varargin{2});
            stoptime = datenum(varargin{3});
            % Take time series from starttime:stoptime
            t_converted = obj.dataset.time(tvar); 
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
          ts.time = obj.dataset.time(tvar); % this is a bad assumption
          % AND---> ALL variables
          vars = obj.dataset.variables;
          for i=1:length(vars)
            dat = obj.dataset.data(vars{i}); % interim notation
            if length(dat)==length(ts.time)
              if ~strcmp(vars{i}, tvar)
                ts.(vars{i}) = [];
                ts.(vars{i}) = dat;
              end
            end
            clear dat
          end
        end
      else % nargin of 3
        % Take time series froms starttime:stoptime
        starttime = datenum(varargin{2});
        stoptime = datenum(varargin{3});
        variables = varargin{3};
        
        t_converted = obj.dataset.time(tvar); % this is a bad assumption
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
    
    function ps = point_size(obj, varName)
      % NCPOINT.point_size - Function to return the size at a point for a specific var.
      % Usage:
      %           size = ncpoint.point_size(varName);
      s = obj.size(varName);
      switch length(s)
        case 1
          ps = s;
        otherwise
          ps = s(2:end);
      end
    end
    
    function c = collection_count(obj, varName)
      % NCPOINT.collection_count - Function to return the number of points in the collection.
      % Usage:
      %           size = ncpoint.collection_count;
      s = obj.size(varName);
      switch length(s)
        case 1
          c = 1;
        otherwise
          c = s(1);
      end
    end
    
    function s = size(obj, varName)
      % NCPOINT.size - Function to return the total size of a specified variable.
      % Usage:
      %           size = ncpoint.size(varName);
      s = obj.dataset.size(varName);
    end
    
    
    
  end % END METHODS
  
  
end % end function