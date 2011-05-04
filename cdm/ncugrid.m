% UGRID Extention of dataset object class for unstructured grid datasets.
% Alexander Crosby, Applied Science Associates 2010
% 
classdef ncugrid < handle
    
  properties (SetAccess = private)
   dataset
   netcdfugrid
   numberofgrids
   meshes
   cells
   nodes
   edges
   faces
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
    
    function ss = unstructuredLatLonSubset(obj, struct) % input subset structure
      % NCUGRID.unstructuredLatLonSubset  - Function to subset an unstructured model grid by lat/lon
      % bounding box using subsetting methods in the ugrid-java.
      ss = []; % Return nothing at the moment, wait for kyle
    end % end subset
  
  end % end methods

end % end class