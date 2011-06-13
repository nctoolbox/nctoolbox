% CGRID Extention of dataset object class for complex/plaid grids (2-D lat,
% lon) often from curvilinear things like ROMS(staggered) or SLOSH. Alexander Crosby,
% Applied Science Associates 2010
% 
% NCTOOLBOX (http://code.google.com/p/nctoolbox)
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
    
  end
  
  
end
