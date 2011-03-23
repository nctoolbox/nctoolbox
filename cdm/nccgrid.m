% CGRID Extention of dataset object class for complex/plaid grids (2-D lat,
% lon) often from curvilinear things like ROMS(staggered) or SLOSH. Alexander Crosby,
% Applied Science Associates 2010
% 
classdef nccgrid < ncsgrid
    
  properties (SetAccess = private)
%     dataset
    
  end
  
  methods
    
    function obj = nccgrid(nc)
      obj = obj@ncsgrid(nc);
%       if ischar(nc)
%         obj.dataset = ncdataset(src);  % nc is a string URL/File
%       elseif isa(nc, 'ncdataset')
%         obj.dataset = nc;             % nc is an ncdataset
%       elseif isa(nc, 'ncsgrid')
%         obj.dataset = nc;            % nc is an ncvariable
%         
%       else
%         ex = MException('POINT:ncobject', 'Incompatable usage, no reference to local netcdf file, opendap url, ncdataset or ncvariable');
%         ex.throw;
%       end
    end
    
    function rg = regrid(obj, lat, lon)
      % Regrid to a specified regular grid
    end
    
    function vr = rotatevec(obj, tvariable)
    end
    
    function gcc = cellcenter(obj)
      % use NJFunc.gridcellcenter(data, grid, type) to regrid both
      % (typically) u/v staggered variables to the grid cell centers, where
      % something like zeta would be located
      
    end
    
    function nt = nearto(obj, nlon, nlat, variable)
      % replicate the point nearto method based on plaid 2d lat/lon grids
      % should return something like indices for just the x,y
      g = obj.dataset.grid(variable);
      v = obj.variable(variable);
      a = v.axes;
      lon = g.(char(a(end)));
      lat = g.(char(a(end-1)));
      
      xdist = abs(lon-nlon); %2d
      ydist = abs(lat-nlat); %2d
      hy_dist = sqrt((xdist.^2)+(ydist.^2)); % 2d matrix of shortest dist
      [row col] = find(hy_dist==min(hy_dist)); 
      obj.where.lon = lon(col);
      obj.where.lat = lat(row);
      obj.where.x = row;
      obj.where.y = col;
      nt = {obj.where.lon obj.where.lat [row col]};
    end
    
    function bb = bbox(obj, t1, t2, z1, z2, east_min, north_min, east_max, north_max, stride)
      vars = obj.dataset.variables;
      for i=1:length(vars)
        A = obj.dataset.variable(vars{i});
        bb.(vars{i}) = A.bbox(t1, t2, z1, z2, east_min, north_min, east_max, north_max, stride);
      end
    end 
    
    function [a b c d] = bboxij(obj, east_min, north_min, east_max, north_max, tvariable)
      % find just the indicies/ranges of grid within the bbox 
      % call dataset.grid(variable)
      % do subset index logic from ncvariable.bbox
      % return those indexs, but how? [lat_min lat_max lon_min lon_max]?
      A = obj.dataset.variable(tvariable);
      [a b c d] = A.bboxij(east_min, north_min, east_max, north_max);
    end
    
  end
  
  
end
