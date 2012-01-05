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
%
% NCTOOLBOX (http://code.google.com/p/nctoolbox)
classdef ncgeovariable < ncvariable
    
    properties (SetAccess = private)
        %         dataset          % ncdataset instance
    end
    
    properties (Dependent = true)
        %         name            % The string variable name that this object represents
        %         attributes\
        %         axes
    end
    
    properties (SetAccess = private, GetAccess = protected)
        %         variable        % ucar.nc2.Variable instance. Represents the data
        %         axesVariables    % ucar.nc2.Variable instance. Represents the data.
        
    end
    
    properties (SetAccess = private, GetAccess = private)
        %         variable        % ucar.nc2.Variable instance. Represents the data
        %         axesVariables    % ucar.nc2.Variable instance. Represents the data.
        axes_info % list of axes names and the dimensions in one cell dict
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
        
        function a = get.axes_info(src)
            switch length(src.size)
                case 1
                    a = fieldnames( src.grid_interop(1,1,1));
                    b = fieldnames( src.grid(1,1,1));
                case 2
                    a = fieldnames( src.grid_interop([1,1],[1,1],[1,1]));
                    b = fieldnames( src.grid([1,1],[1,1],[1,1]));
                case 3
                    a = fieldnames( src.grid_interop([1,1,1], [1,1,1], [1,1,1]));
                    b = fieldnames( src.grid([1,1,1], [1,1,1], [1,1,1]));
                case 4
                    a = fieldnames( src.grid_interop([1,1,1,1], [1,1,1,1], [1,1,1,1]));
                    b = fieldnames( src.grid([1,1,1,1], [1,1,1,1], [1,1,1,1]));
                case 5
                    a = fieldnames( src.grid_interop([1,1,1,1,1], [1,1,1,1,1], [1,1,1,1,1]));
                    b = fieldnames( src.grid([1,1,1,1,1], [1,1,1,1,1], [1,1,1,1,1]));
                case 6
                    a = fieldnames( src.grid_interop([1,1,1,1,1,1], [1,1,1,1,1,1], [1,1,1,1,1,1]));
                    b = fieldnames( src.grid([1,1,1,1,1,1], [1,1,1,1,1,1], [1,1,1,1,1,1]));
            end
            len = length(a);
            for i = 1:len
                a{i, 2} = src.dataset.size(b{i});
            end
            
        end
        
        function e = extent(src)
            % NCGEOVARIABLE.extent - Function to find geographic bounding box coordinates.
            % Usage e = var.extent;
             s = src.size;
             lens = length(s);
             first = ones(1,lens);
             if lens > 1
                 s(1) = 1;
                 if lens > 3;
                     s(2) = 1;
                 end
             end
             stride = first;
%             switch length(s) % hopefully this speeds up the grid_interop call when time is involved
%                 case 1
%                     g = src.grid_interop(:);
%                 case 2
%                     g = src.grid_interop(:,:);
%                 case 3
%                     g = src.grid_interop(1,:,:);
%                 case 4
                    g = src.grid_interop(first, s, stride);
%             end
            e.lon = [min(min(g.lon)) max(max(g.lon))];
            e.lat = [min(min(g.lat)) max(max(g.lat))];
        end % end extent
        
        function te = timeextent(src)
            % NCGEOVARIABLE.timeextent - Function to find the start and stop time of the variable.
            % Usage e = var.timeextent;
            s = src.size;
            lens = length(s);
            first = ones(1, lens);
            stride = first;
            last = first;
            last(1) = s(1);
            
            %                     g1 = src.grid_interop(first, first, stride);
            %                     g2 = src.grid_interop(s, s, stride);
            g = src.grid_interop(first, last, stride);
            
            %             te = [g1.time g2.time];
            if isfield(g, 'time')
                te = [min(g.time) max(g.time)];
            else
                error('NCGEOVARIABLE:TIMEEXTENT',...
                    'There appears to be no time axis associated with the variable.');
            end
        
        end % end timeextent
        
        function ig = grid_interop(src, first, last, stride)
            % NCGEOVARIABLE.GRID_INTEROP - Method to get the coordinate variables and their data as a
            % a structure with standardized field names for lat, lon, time, and z. Other coordiante variables
            % that are not recognized by netcdf-java as the previous four types have field names directly
            % taken from their variable names.
            % Useage: >> gridstruct = geovar.grid_interop(1,:,:,1:2:50);
            
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
                            if strcmpi(pos_z, 'POSITIVE_DOWN')
                                tmp = g.(tempname);
                                ig.z = tmp.*-1; %adjust for positive direction
                            elseif strcmpi(pos_z, 'down')
                                tmp = g.(tempname);
                                ig.z = tmp.*-1; %adjust for positive direction
                            else
                                ig.z = g.(tempname);
                            end
                            
                        case 'GeoZ'
                            pos_z = char(javaaxisvar.getPositive());
                            z_sn = src.dataset.attribute(tempname, 'standard_name');
                            k = strfind(z_sn, 'ocean_s');
                            switch isempty(k)
                                case 0
                                    %                       n = max(size(src.size));
                                    %                       trange = java.util.ArrayList(n);
                                    %                       zrange = trange;
                                    %                       xrange = trange;
                                    %                       yrange = trange;
                                    trange = ucar.ma2.Range(first(1) - 1, last(1) - 1, stride(1));
                                    zrange = ucar.ma2.Range(first(2) - 1, last(2) - 1, stride(2));
                                    xrange = ucar.ma2.Range(first(4) - 1, last(4) - 1, stride(4));
                                    yrange = ucar.ma2.Range(first(3) - 1, last(3) - 1, stride(3));
                                    %                       coordinates.GridCoordSys.getVerticalTransform.getCoordinateArray
                                    %                          temp = src.variable.getCoordinateSystems();
                                    griddataset = ucar.nc2.dt.grid.GridDataset.open(src.dataset.location);
                                    grid = griddataset.findGridByName(src.name);
                                    grid = grid.getCoordinateSystem();
                                    subgrid = grid.getVerticalTransform();
                                    subgrid = subgrid.subset(trange, zrange, yrange, xrange); %tempfortesting arc 10/18
                                    %
                                    % It looks like this dataset (http://geoport.whoi.edu/thredds/dodsC/usgs/vault0/models/examples/bora_feb.nc)
                                    % works when the vertical transform isnt subset. NJTBX uses the same methodology but subsets the resulting
                                    % z coordinate field in matlab after calling the getCoordinateArray method, not before in the java.
                                    % Should I implement a try and then default to subsetting afterwards on catch? Or should I just always do it
                                    % after? I am choosing the former. Anyone have any thoughts on the matter? -acrosby
                                    try
                                        try
                                            try
                                                for q = first(1):stride(1):last(1)
%                                                     grid = griddataset.findGridByName(src.name);
%                                                     subgrid = grid.getCoordinateSystem();
                                                    array = subgrid.getCoordinateArray(q-1);
                                                    ig.z(q, :, :, :) = array.copyToNDJavaArray();
                                                end
                                            catch me
                                                c = 1;
                                                for q = first(1):stride(1):last(1)
%                                                     grid = griddataset.findGridByName(src.name);
%                                                     grid = grid.getCoordinateSystem();
                                                    subgrid = grid.getVerticalTransform();
                                                    
                                                    array = subgrid.getCoordinateArray(q-1); % Issue 27 is failing here...
                                                    z(c, :, :, :) = array.copyToNDJavaArray();
                                                    c = c + 1;
                                                end
                                                ig.z = z(:, first(2):stride(2):last(2),  first(3):stride(3):last(3),  first(4):stride(4):last(4));
                                            end
                                        catch me
                                            array = subgrid.getCoordinateArray(0);
                                            ig.z = array.copyToNDJavaArray();
                                        end
                                    catch me                                        
                                        disp('Could you please add the code you are trying to run to Issue 27 at the nctoolbox issue tracking site.');
                                        web http://code.google.com/p/nctoolbox/issues/detail?id=27
                                        me.throw()
                                        % me.error('There is a problem applying the vertical coordinate tranform and subsetting the resuting values.');
                                    end
                                    
                                otherwise
                                    
                                    if strcmp(pos_z, 'POSITIVE_DOWN')
                                        tmp = g.(tempname);
                                        ig.z = tmp.*-1; %adjust for positive direction
                                    else
                                        ig.z = g.(tempname);
                                    end
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
%                             ind = find(tmp > 180); % convert 0-360 convention to -180-180
%                             tmp(ind) = tmp(ind)-360;
                            ig.lon = tmp;    
                        case 'Lat'
                            ig.lat = g.(tempname);
                            
                        case 'GeoY'
                            ig.y = g.(tempname);
                            if exist('griddataset', 'var')
                                grid = griddataset.findGridByName(src.name);
                                grid = grid.getCoordinateSystem();
                            else
                                griddataset = ucar.nc2.dt.grid.GridDataset.open(src.dataset.location);
                                grid = griddataset.findGridByName(src.name);
                                grid = grid.getCoordinateSystem();
                            end
                            try
                                %ig.y = grid.getYHorizAxis;
                                %ig.x = grid.getXHorizAxis;
                                [x, y] = meshgrid(g.x,g.y);
                                s = size(x);
                                x = reshape(x, [1 numel(x)]);
                                y = reshape(y, [1 numel(y)]);
                                tempXY = [x; y];
                                projection = grid.getProjection();
                                tempLatLon = projection.projToLatLon(tempXY);
                                ig.lat = reshape(tempLatLon(1,:), s); 
                                ig.lon = reshape(tempLatLon(2,:), s); 
                                
                            catch me
                            end

%                         case 'GeoX'
%                             ig.x = g.(tempname);
%                             if exist('griddataset', 'var')
%                                 grid = griddataset.findGridByName(src.name);
%                                 grid = grid.getCoordinateSystem();
%                             else
%                                 griddataset = ucar.nc2.dt.grid.GridDataset.open(src.dataset.location);
%                                 grid = griddataset.findGridByName(src.name);
%                                 grid = grid.getCoordinateSystem();
%                                 
%                             end
%                             try
%                                 ig.x = grid.getXHorizAxis;
%                                 projection = grid.getProjection();
%                                 ig.lon = projection.projToLatLon();
%                             catch me
%                             end

                
                        otherwise
                            ig.(tempname) = g.(tempname);
                            
                    end % end switch on type
                end % end is type empty or not if statement
            end % end loop through field names
            
        end % grid_interop end
        
        function tw = timewindow(src, varargin)
            % NCGEOVARIABLE.TIMEWINDOW - Function to pull the time coordinates within the specified
            % start and stop times from the variable object.
            % Useage: >> time = geovar.timewindow([2004 1 1 0 0 0], [2005 12 31 0 0 0]);
            %              >> time = geovar.timewindow(731947, 732677);
            if nargin < 3
                d = src.timewindowij(varargin{1});
            elseif nargin == 3
                d = src.timewindowij(varargin{1}, varargin{2});
            else
                error('NCGEOVARIABLE:TIMEWINDOW',...
                    'Too many input arugments.');
            end
            tw = d;
        end
        
        function tv = gettimevar(src)
            % NCGEOVARIABLE.gettimevar()
            tn = src.gettimename();
            tv = src.dataset.geovariable(tn);
        end
        
        function lv = getlonvar(src)
            % NCGEOVARIABLE.getlonvar()
            tn = src.getlonname();
            lv = src.dataset.geovariable(tn);
        end
        
        function lv = getlatvar(src)
            % NCGEOVARIABLE.getlatvar()
            tn = src.getlatname();
            lv = src.dataset.geovariable(tn);
        end
        
        function tn = gettimename(src)
            % NCGEOVARIABLE.gettimename()
            for i = 1:length(src.axes)
                tempname = src.axes{i};
                javaaxisvar  =   src.dataset.netcdf.findVariable(tempname);
                type{i} = char(javaaxisvar.getAxisType());
            end
            match = strcmp('Time', type);
            tn = src.axes(match);
        end
        
        function ln = getlonname(src)
            % NCGEOVARIABLE.getlonname()
            for i = 1:length(src.axes)
                tempname = src.axes{i};
                javaaxisvar  =   src.dataset.netcdf.findVariable(tempname);
                type{i} = char(javaaxisvar.getAxisType());
            end
            match = strcmp('Lon', type);
            ln = src.axes(match);
        end
        
        function ln = getlatname(src)
            % NCGEOVARIABLE.gelatname()
            for i = 1:length(src.axes)
                tempname = src.axes{i};
                javaaxisvar  =   src.dataset.netcdf.findVariable(tempname);
                type{i} = char(javaaxisvar.getAxisType());
            end
            match = strcmp('Lat', type);
            ln = src.axes(match);
        end
        
         function tn = gettimedata(src, start, last, stride)
            % NCGEOVARIABLE.gettimedata()
            var = src.gettimevar;
            tn = var.data(start, last, stride);
            tn = var.dataset.time(src.gettimename, tn);
        end
        
        function s = getlondata(src, start, last, stride)
            % NCGEOVARIABLE.getlondata()
            v = src.getlonvar;
            sz = src.size();
            lonsize = v.size();
            lonlocation = find(sz==lonsize(1));
             if length(lonsize) == 1
                lonstart = start(lonlocation);
                lonlast = last(lonlocation);
                lonstride = stride(lonlocation);
            else
                lonstart = start(lonlocation:lonlocation+1);
                lonlast = last(lonlocation:lonlocation+1);
                lonstride = stride(lonlocation:lonlocation+1);
            end
            switch length(lonsize)
                  case 1
                    s = v.data(lonstart:lonstride:lonlast);
                  case 2
                    s = v.data(lonstart(1):lonstride(1):lonlast(1),lonstart(2):lonstride(2):lonlast(2));
                  case 3
                    s = v.data(lonstart(1):lonstride(1):lonlast(1), lonstart(2):lonstride(2):lonlast(2), lonstart(3):lonstride(3):lonlast(3));
                  case 4
                    s = v.data(lonstart(1):lonstride(1):lonlast(1), lonstart(2):lonstride(2):lonlast(2), lonstart(3):lonstride(3):lonlast(3),...
                      lonstart(4):lonstride(4):lonlast(4));
            end
        end
        
        function s = getlatdata(src, start, last, stride)
            % NCGEOVARIABLE.gelatdata()
            v = src.getlatvar;
            sz = src.size();
            latsize = v.size();
            latlocation = find(sz==latsize(1));
            if length(latsize) == 1
                latstart = start(latlocation);
                latlast = last(latlocation);
                latstride = stride(latlocation);
            else
                latstart = start(latlocation:latlocation+1);
                latlast = last(latlocation:latlocation+1);
                latstride = stride(latlocation:latlocation+1);
            end
            switch length(latsize)
                  case 1
                    s = v.data(latstart:latstride:latlast);
                  case 2
                    s = v.data(latstart(1):latstride(1):latlast(1),latstart(2):latstride(2):latlast(2));
                  case 3
                    s = v.data(latstart(1):latstride(1):latlast(1), latstart(2):latstride(2):latlast(2), latstart(3):latstride(3):latlast(3));
                  case 4
                    s = v.data(latstart(1):latstride(1):latlast(1), latstart(2):latstride(2):latlast(2), latstart(3):latstride(3):latlast(3),...
                      latstart(4):latstride(4):latlast(4));
            end
        end
        
        
        
        %% These functions would rather output multiple outputs instead of struct, must reconcile
        %     with the subsref in either ncgeovariable or ncvariable. Wait, why, then, does geoij work???
        function d = timewindowij(src, varargin)
            % NCGEOVARIABLE.TIMEWINDOWIJ - Function to get indices from start and stop times for sub-
            % setting. TODO: There must be a better/fast way to do this using the java library.
            % Useage: >> timestruct = geovar.timewindowij([2004 1 1 0 0 0], [2005 12 31 0 0 0]);
            %              >> timestruct = geovar.timewindowij(731947, 732677);
            % Result: time.time = time data
            %            time.index = time indices
            s = src.size;
            first = ones(1, length(s));
            last = s;
            stride = first;
            g = src.grid_interop(first, last, stride);
            
            if isfield(g, 'time') % are any of the fields recognized as time explictly
                if nargin > 2 % If two times are input, do window
                    starttime = datenum(varargin{1});
                    stoptime = datenum(varargin{2});
                    if isempty(starttime)
                        starttime = g.time(1);
                    end
                    if isempty(stoptime)
                        stoptime = g.time(end);
                    end
                    
                    t_index1 = g.time >= starttime;
                    t_index2 = g.time <= stoptime;
                    d.index = find(t_index1==t_index2);
                    d.time = g.time(d.index);
                elseif nargin < 3 % If theres only one time, do nearest
                    neartime = datenum(varargin{1});
                    diff = abs(g.time-neartime);
                    ind = find(diff==min(diff));
                    d.index = ind(1);
                    d.time = g.time(d.index);
                    if length(ind) > 1;
                        warning('NCGEOVARIABLE:TIMEWINDOWIJ',...
                            ['Multiple time indices determined to be nearest to the supplied time,'...
                            'only the first index was output.']);
                    end
                end
            else
                me = MException(['NCTOOLBOX:ncgeovariable:timewindowij'], ...
                    'No grid variable returned as time.');
                me.throw;
            end
        end % end timewindowij
        
        function d = geosubset(obj, struct)
            % NCGEOVARIABLE.GEOSUBSET -
            % Create subset structure:
            % s.time=[now-3 now]; % Can be omitted or s.t_index=[i1 i2]; can be used for subsetting time using indices
            % s.t_stride=2; % Can be omitted
            % s.lat=[22 27];
            % s.lon=[117 123];
            % s.h_stride=[2 2]; % Can be omitted
            % s.z_index=[1 30]; % Can be omitted
            % s.v_stride=2; % Can be omitted
            %
            % Usage:
            % sub = gvar.geosubset(s); % Subset method
            %
            % Returns:
            % Structure with fields 'data' for data values, and 'grid' for structure of coordinate variables
            % (that come from grid_interop).
            
            
            %           if ~regexp(obj.dataset.attribute('Conventions'), 'UGRID')
            nums = obj.size;
            if isfield(struct, 'h_stride');
                if length(struct.h_stride) < 2
                    struct.h_stride(1, 2) = struct.h_stride(1,1); % Resoltution to Issue 24
                end
            else
                struct.h_stride = [1 1];
            end
            
            if isfield(struct, 'v_stride');
            else
                struct.v_stride = 1;
            end
            
            if isfield(struct, 't_stride');
            else
                struct.t_stride = 1;
            end
            
            [indstart_r indend_r indstart_c indend_c] = obj.geoij(struct);
            
            if ~isempty(indstart_r)
                if isfield(struct, 'time') % Deal with time (values) or t_index (indices) bounds
                    if iscell(struct.time)
                        switch length(struct.time)
                            case 1
                                t = obj.timewindowij(struct.time{1});
                                tmin_i = t.index;
                                tmax_i = t.index;
                            case 2
                                t = obj.timewindowij(struct.time{1}, struct.time{2});
                                tmin_i = min(t.index);
                                tmax_i = max(t.index);
                        end
                    else
                        if ischar(struct.time)
                            t = obj.timewindowij(struct.time);
                            tmin_i = t.index;
                            tmax_i = t.index;
                        elseif isarray(struct.time)
                            t = obj.timewindowij(struct.time);
                            tmin_i = t.index;
                            tmax_i = t.index;
                        else
                            switch length(struct.time)
                                case 1
                                    t = obj.timewindowij(struct.time);
                                    tmin_i = t.index;
                                    tmax_i = t.index;
                                case 2
                                    t = obj.timewindowij(struct.time(1), struct.time(2));
                                    tmin_i = min(t.index);
                                    tmax_i = max(t.index);
                                otherwise % for anything else assume that it is a single datevec
                                    t = obj.timewindowij(struct.time);
                                    tmin_i = t.index;
                                    tmax_i = t.index;
                            end
                        end
                    end
                elseif isfield(struct, 't_index') % check for 1 time index not the same time index twice $$$$$$$$$$$$
                    if iscell(struct.t_index)
                        switch length(t_index)
                            case 1
                                if numel(struct.t_index{1}) > 1 % check to see if someone used str or datevec by accident
                                    me = MException(['NCTOOLBOX:ncgeovariable:geosubset'], ...
                                        'Expected min time to be an index/integer.');
                                    me.throw;
                                else
                                    tmin_i = struct.t_index{1};
                                end
                            case 2
                                if numel(struct.t_index{1}) > 1 % check to see if someone used str or datevec by accident
                                    me = MException(['NCTOOLBOX:ncgeovariable:geosubset'], ...
                                        'Expected min time to be an index/integer.');
                                    me.throw;
                                else
                                    tmin_i = struct.t_index{1};
                                end
                                if numel(struct.t_index{2}) > 1 % check to see if someone used str or datevec by accident
                                    me = MException(['NCTOOLBOX:' mfilename ':geosubset'], ...
                                        'Expected max time to be an index/integer.');
                                    me.throw;
                                else
                                    tmax_i = struct.t_index{2};
                                end
                        end
                    else
                        switch length(struct.t_index)
                            case 1
                                tmin_i = struct.t_index(1);
                                tmax_i = struct.t_index(1);
                            case 2
                                tmin_i = struct.t_index(1);
                                tmax_i = struct.t_index(2);
                        end
                    end
                else
                    tmin_i = 1;
                    tmax_i = nums(1);
                end
                
                order = obj.getaxesorder;
                
                if isfield(struct, 'z_index')
                    switch length(struct.z_index)
                        case 1
                            zmin = struct.z_index;
                            zmax = struct.z_index;
                        case 2
                            zmin = struct.z_index(1);
                            zmax = struct.z_index(2);
                    end
                else
                    if isfield(order, 'z')
                        zmin = 1;
                        zmax = nums(order.z);
                    else
                        order.z = [];
                        zmin = [];
                        zmax = [];
                    end
                end

                if length(nums) < 2
                    me = MException(['NCTOOLBOX:ncgeovariable:geosubset'], ...
                        ['Expected data of ', obj.name, ' to be at least rank 2.']);
                    me.throw;
                else
                    first = ones([1, length(nums)]);
                    stride = ones([1, length(nums)]);
                    last = nums;
                    
                    if order.lon ~= order.lat      
                        stride(order.time)   = struct.t_stride; 
                        stride(order.z)         = struct.v_stride;
                        stride(order.lon)     = struct.h_stride(2);
                        stride(order.lat)      = struct.h_stride(1);
                    else
                        order.lat = order.lon + 1;
                        stride(order.time)   = struct.t_stride;
                        stride(order.z)         = struct.v_stride;
                        stride(order.lon)     = struct.h_stride(1);
                        stride(order.lat)      = struct.h_stride(2);
                    end
                    first(order.time)   = tmin_i;
                    first(order.z)         = zmin;
                    first(order.lon)     = indstart_r;
                    first(order.lat)      = indstart_c;
                    last(order.time)   = tmax_i;
                    last(order.z)         = zmax;
                    last(order.lon)     = indend_r;
                    last(order.lat)      = indend_c;
                    
                    %                 else
                    %                     me = MException(['NCTOOLBOX:ncgeovariable:geosubset'], ...
                    %                         ['Expected shape of data of ', obj.name, ' to be less than rank 5.']);
                    %                     me.throw;
                end
                
                % Get the corresponding data and interop grid...
                d.data = obj.data(first, last, stride);
                d.grid = obj.grid_interop(first, last, stride);
                %           else
                %             ugrid = ncugrid(obj); % Starting to add place holders for ugrid subsetting functionality
                %             d = ugrid.unstructuredLatLonSubset(struct);
                %           end % end of ugrid if
            else
                me = MException(['NCTOOLBOX:ncgeovariable:geosubset'], ...
                        ['No indices returned cooresponding to geographic subset.']);
                me.throw;
            end
        end % end of geosubset
        %%
        function [indstart_r indend_r indstart_c indend_c] =...
                geoij(obj, struct)
            % GEOVARIABLE.GEOIJ - Function to return start and stop indices for rows and columns if the
            % the lat/lon grids are rank 2, but if the grids are vectors the function returns start and stop
            % indices of lon and then lat (in that order).
            %
            % Usage: >> [firstrow, lastrow, firstcol, lastcol] = geoij(geovar, subsetstruct)
            %            >> [firstlon, lastlon, firstlat, lastlat] = geoij(geovar, subsetstruct) % if lat/lon are vector
            %
            %
            %             s = obj.size;
            %             first = ones(1, length(s));
            %             last = s;
            %             stride = first;
            %             g = obj.grid_interop(first, last, stride); % this is doing the whole thing
            %           h = 0;
            flag = 0;
            
            warning off
            switch  length(obj.size)
                case 1
                    g.lon = obj.getlondata(1, obj.size, 1);
                    g.lat = obj.getlatdata(1, obj.size, 1);
                case 2
                    g.lon = obj.getlondata([1, 1], obj.size, [1, 1]);
                    g.lat = obj.getlatdata([1, 1], obj.size, [1, 1]);
                case 3
                    g.lon = obj.getlondata([1, 1, 1], obj.size, [1, 1, 1]);
                    g.lat = obj.getlatdata([1, 1, 1], obj.size, [1, 1, 1]);
                case 4
                    g.lon = obj.getlondata([1, 1, 1, 1], obj.size, [1, 1, 1, 1]);
                    g.lat = obj.getlatdata([1, 1, 1, 1], obj.size, [1, 1, 1, 1]);
            end
            warning on
            
            if max(g.lon) > 360
                if min(g.lon) <~ 0
                    ind = find(g.lon < 0); % convert 0-360 convention to -180/180
                    g.lon(ind) = g.lon(ind)+360;
                else
                    error('NCGEOVARIABLE:GEOIJ',...
                        'Longitude contains values that follow both -180/180 and 0/360+ conventions; can not subset.');
                end
            elseif min(g.lon) <~ 0
                % Do nothing, we assume that input is in -180/180 conventions. And this means that the data is too.
            end
            
            
            %Unpack geosubset_structure
            if isfield(struct, 'lat');
                switch length(struct.lat)
                    case 1
                        flag = 1;
                    case 2
                        north_max = struct.lat(2);
                        north_min = struct.lat(1);
                end
                
            else
                north_max = max(g.lat);
                north_min = min(g.lat);
            end
            
            if isfield(struct,'lon');
                switch length(struct.lon)
                    case 1
                        flag = 1;
                    case 2
                        east_max = struct.lon(2);
                        east_min = struct.lon(1);
                end
                
            else
                east_max = max(g.lon);
                east_min = min(g.lon);
            end
            
            switch flag
                case 0
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
                        indstart_c = min(ind_c); % Out
                        indend_c = max(ind_c);
                        indstart_r = min(ind_r);
                        indend_r = max(ind_r);
                    else
                        indlat1 = (g.lat <= north_max);
                        indlat2 = (g.lat >= north_min);
                        indlat = find(indlat1&indlat2);
                        indlon1 = (g.lon <= east_max);
                        indlon2 = (g.lon >= east_min);
                        indlon = find(indlon1&indlon2);
                        % Out
                        attrs = obj.dataset.attributes;
                        key = value4key(attrs, 'CF:featureType');
                        
                        if strcmp(key, 'timeSeries')
                            a = [];
                            for i = 1:length(indlon)
                                if ~isempty(find(indlon(i)==indlat, 1))
                                    a(length(a)+1) = indlon(i);
                                end
                            end
                            indstart_c = min(a);
                            indend_c = max(a);
                            indstart_r = min(a);
                            indend_r = max(a);
                        else
                            indstart_c = min(indlon); % testing temp, i switched lon and lat here
                            indend_c = max(indlon);
                            indstart_r = min(indlat);
                            indend_r = max(indlat);
                            %
                        end
                    end
                case 1
                    if isvector(g.lat)
                        indstart_c = near(g.lon,struct.lon);
                        indstart_r = near(g.lat,struct.lat);
                        indend_c = indstart_c;
                        indend_r = indstart_r;
                    else
                        % This way does't require the addtional step, or the additional function and is verified
                        % to return the same results as using nearxy with first output and ind2ij.
                        % Ind2ij doesn't seem to add any value to the toolbox, everything is already here to
                        % get rows and columns. The third nearxy output can remain private if it confusing to
                        % users (no mention of it in help docs).
                        [a b indexes] = nearxy(g.lon, g.lat, struct.lon, struct.lat);
                        indstart_c = indexes(2);
                        indend_c = indexes(2);
                        indstart_r = indexes(1);
                        indend_r = indexes(1);
                    end
            end
        end % end geoij
        
        function order = getaxesorder(obj)
            %             permute_nums = fliplr(obj.axes_info);
            ainfo = obj.axes_info;
            siz = obj.size;
            for i = 1:length(ainfo(:, 1))
                order.(ainfo{i,1}) = find(siz==ainfo{i, 2}(1));
            end
            
        end
        
        function sref = subsref(obj,s)
            switch s(1).type
                % Use the built-in subsref for dot notation
                case '.'
                    switch s(1).subs
                        case 'gettimedata'
                            switch length(s)
                                case 1
                                    sref = obj;
                                case 2
                                    nums = obj.size;
                                    [first last stride] = indexing(s(2).subs, double(nums));
                                    sref = obj.gettimedata(first, last, stride);
                            end
                            case 'getlondata'
                            switch length(s)
                                case 1
                                    sref = obj;
                                case 2
                                    nums = obj.size;
                                    [first last stride] = indexing(s(2).subs, double(nums));
                                    sref = obj.getlondata(first, last, stride);
                            end
                            case 'getlatdata'
                            switch length(s)
                                case 1
                                    sref = obj;
                                case 2
                                    nums = obj.size;
                                    [first last stride] = indexing(s(2).subs, double(nums));
                                    sref = obj.getlatdata(first, last, stride);
                            end
                        case 'grid_interop'
                            switch length(s)
                                case 1
                                    sref = obj;
                                case 2
                                    nums = obj.size;
                                    [first last stride] = indexing(s(2).subs, double(nums));
                                    sref = obj.grid_interop(first, last, stride);
                            end
                        case 'data'
                            nums = size(obj);
                            if ~isempty(nums)
                                switch length(s)
                                    case 1
                                        sref = obj;
                                    case 2
                                        [first last stride] = indexing(s(2).subs, double(nums));
                                        sref = obj.data(first, last, stride);
                                end
                                
                            else
                                sref = obj.data;
                                warning(['NCTOOLBOX:ncgeovariable:subsref'], ...
                                    ['Variable "' name '" has no netcdf dimension associated with it. Errors may result from non CF compliant files.'])
                            end
                        case 'grid'
                            nums = size(obj);
                            if ~isempty(nums)
                                switch length(s)
                                    case 1
                                        sref = obj;
                                    case 2
                                        [first last stride] = indexing(s(2).subs, double(nums));
                                        sref = obj.grid(first, last, stride);
                                end
                                
                            else
                                warning(['NCTOOLBOX:ncgeovariable:subsref'], ...
                                    ['Variable "' name '" has no netcdf dimension associated with it. Errors may result from non CF compliant files.'])
                            end
                        otherwise
                            sref = builtin('subsref',obj,s);
                    end
                case '()'
                    warning(['NCTOOLBOX:ncgeovariable:subsref'], ...
                        'Not a supported subscripted reference, "()" are not permitted to call variable object methods');
                case '{}'
                    warning(['NCTOOLBOX:ncgeovariable:subsref'], ...
                        'Not a supported subscripted reference, "{}" are not permitted to call variable object methods');
            end
        end
    end % methods end
    
    %     methods (Access = protected)
    %
    %     end % protected methods end
end % class end