% SGRID Extention of dataset object class for structured grids.
% Alexander Crosby, Applied Science Associates 2010
% 
classdef ncsgrid < handle
    
  properties (SetAccess = private)
    dataset        
  end
  
  methods
    
    function obj = ncsgrid(nc)
      if ischar(nc)
        obj.dataset = ncdataset(src);  % nc is a string URL/File
      elseif isa(nc, 'ncdataset')
        obj.dataset = nc;             % nc is an ncdataset
%       elseif isa(nc, 'ncsgrid')
%         obj.dataset = nc.dataset;            % nc is an ncvariable
        
      else
        ex = MException('SGRID:ncobject', 'Incompatable usage, no reference to local netcdf file, opendap url, ncdataset or ncvariable');
        ex.throw;
      end

    end
    
    function gt = gridtype(obj, variable)
      g = obj.dataset.grid(variable);
      lat = g.lat;
      
      if isvector(lat);
        gt = 'rgrid';
      else
        gt = 'cgrid';
      end
      
    end

    function rg = regrid(obj, variable, lat, lon)
      % SGRID.regrid = Function to regrid data of a dataset object to an input grid.
    rg = []
    end
    
    function r = rgrid(obj)
      r = ncrgrid(obj);
    end
    
    function c = cgrid(obj)
      c = nccgrid(obj)
    end
    
    
  end
  
  
end