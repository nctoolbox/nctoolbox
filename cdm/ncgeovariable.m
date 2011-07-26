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
                            z_sn = src.dataset.attribute('standard_name', tempname);
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
                                    grid = ucar.nc2.dt.grid.GridDataset.open(src.dataset.location);
                                    grid = grid.findGridByName(src.name);
                                    grid = grid.getCoordinateSystem();
                                    subgrid = grid.getVerticalTransform();
                                    subgrid = subgrid.subset(trange, zrange, yrange, xrange);
                                    %                            try
                                    array = subgrid.getCoordinateArray(0);
                                    %                            catch
                                    %                              array = subgrid.getCoordinateArray();
                                    %                            end
                                    
                                    ig.z = array.copyToNDJavaArray();
                                    if strcmp(pos_z, 'POSITIVE_DOWN')
                                        tmp = g.z;
                                        ig.z = tmp.*-1; %adjust for positive direction
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
            
            if isfield(struct, 'time') % Deal with time (values) or t_index (indices) bounds
                if iscell(struct.time)
                    t = obj.timewindowij(struct.time{1}, struct.time{2});
                    tmin_i = min(t.index);
                    tmax_i = max(t.index);
                else
                    t = obj.timewindowij(struct.time(1), struct.time(2));
                    tmin_i = min(t.index);
                    tmax_i = max(t.index);
                end
            elseif isfield(struct, 't_index')
                if iscell(struct.t_index)
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
                else
                    tmin_i = struct.t_index(1);
                    tmax_i = struct.t_index(2);
                end
            else
                tmin_i = 1;
                tmax_i = nums(1);
            end
            
            
            
            if length(nums) < 2
                me = MException(['NCTOOLBOX:ncgeovariable:geosubset'], ...
                    ['Expected data of ', obj.name, ' to be at least rank 2.']);
                me.throw;
            elseif length(nums) < 3
                ax = obj.grid([1 1],[1 1],[1 1]);
                if isfield(ax, 'time')
                    first = [tmin_i indstart_r];
                    last = [tmax_i indend_r];
                    stride = [struct.t_stride struct.h_stride(2)];
                else
                    first = [indstart_r indstart_c];
                    last = [indend_r indend_c];
                    stride = [struct.h_stride(2) struct.h_stride(1)];
                end
            elseif length(nums) < 4
                ax = obj.grid([1 1 1],[1 1 1],[1 1 1]);
                if isfield(ax, 'time')
                    first = [tmin_i indstart_r indstart_c];
                    last = [tmax_i indend_r indend_c];
                    stride = [struct.t_stride struct.h_stride(2) struct.h_stride(1)];
                elseif isfield(ax, 'z')
                    if isfield(struct, 'z_index');
                    else
                        struct.z_index = [1 nums(1)];
                    end
                    first = [struct.z_index(1) indstart_r indstart_c];
                    last = [struct.z_index(2) indend_r indend_c];
                    stride = [struct.v_stride struct.h_stride(2) struct.h_stride(1)];
                else
                    me = MException(['NCTOOLBOX:ncgeovariable:geosubset'], ...
                        'Expected either a coordinate variable acknowleged as time or as z.');
                    me.throw;
                end
            elseif length(nums) < 5
                if isfield(struct, 'z_index');
                else
                    struct.z_index = [1 nums(2)];
                end
                first = [tmin_i struct.z_index(1) indstart_r indstart_c];
                last = [tmax_i struct.z_index(2) indend_r indend_c];
                stride = [struct.t_stride struct.v_stride struct.h_stride(2) struct.h_stride(1)];
            else
                me = MException(['NCTOOLBOX:ncgeovariable:geosubset'], ...
                    ['Expected data of ', obj.name, ' to be less than rank 5.']);
                me.throw;
                
            end
            
            % Get the corresponding data and interop grid...
            d.data = obj.data(first, last, stride);
            d.grid = obj.grid_interop(first, last, stride);
            %           else
            %             ugrid = ncugrid(obj); % Starting to add place holders for ugrid subsetting functionality
            %             d = ugrid.unstructuredLatLonSubset(struct);
            %           end % end of ugrid if
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
            s = obj.size;
            first = ones(1, length(s));
            last = s;
            stride = first;
            g = obj.grid_interop(first, last, stride);
            %           h = 0;
            flag = 0;
            
            
            %Unpack geosubset_structure
            if isfield(struct, 'lat');
                north_max = struct.lat(2);
                north_min = struct.lat(1);
            else
                north_max = max(g.lat);
                north_min = min(g.lat);
            end
            
            if isfield(struct,'lon');
                east_max = struct.lon(2);
                east_min = struct.lon(1);
            else
                east_max = max(g.lon);
                east_min = min(g.lon);
            end
            
            
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
                    flag = 1;
                end
            end
            
            
            
        end % end geoij
        
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