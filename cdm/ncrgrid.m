% RGRID Extention of dataset object class for regular grids (1-D lat, lon).
% Alexander Crosby, Applied Science Associates 2010
% 
classdef ncrgrid < ncsgrid
    
  properties (SetAccess = private)
%     dataset              
  end
  
  methods
    function obj = ncrgrid(nc)
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
    
    function nt = nearto(obj, nlat, nlon, variable)
      % replicate the point nearto method based on regular grids 1d lat/lon
      % should return something like indices for just the x,y
      g = obj.dataset.grid(variable);
      
      [latname lonname] = cnames(variable)
%       if isfield(g, 'lat_u')
%           lat = g.lat_u;
%           lon = g.lon_u;
%         
%       elseif isfield(g, 'Y')
%           lat = g.Y;
%           lon = g.X;
%         
%       elseif isfield(g, 'y')
%           lat = g.y;
%           lon = g.x;
       
%       else
%         error('bbox:dimensions','No dimension handle found, specifically latitude')
%       end
      
      xdist = abs(lon-nlon);
      ydist = abs(lat-nlat);
      [d2_xdist d2_ydist] = meshgrid(xdist, ydist);
      hy_dist = sqrt((d2_xdist.^2)+(d2_ydist.^2)); % should be a 2d shortest distance matrix
      mindistmat = ones(size(hy_dist)).*min(min(hy_dist));
      [row col] = find(hy_dist==mindistmat); 
%       obj.where.lon = lon(col);
%       obj.where.lat = lat(row);
%       obj.where.x = row;
%       obj.where.y = col;
      nt = {lat(row) lon(col) [row col]};
    end % fix this
    
    
    
  end
  
  
end
