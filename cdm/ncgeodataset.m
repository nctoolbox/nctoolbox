% NCGEODATASET  Provide access to CF/COARDS convention datasets accessable by the
% NetCDF 4 API. CFDATASET is a subclass of NCDATASET
%
% Use as:
%   ds = ncgeodataset(dataref)
%
% Arguments:
%   dataref = A reference to a ncdataset that can be accessed by the NetCDF 4
%       API. This includes local netcdf files, netcdf files on web servers
%       and OpenDAP URLs
%
% Return:
%   An instance of a ncgeodataset class
%   !!! See ncgeodataset and cfdataset for additional information !!!
%
% Methods:
%   ncgeodataset.grid - Retrieve coordinate data for the variable.
%   ncgeodataset.struct - Retrieve data and coordinate data for the variable
%   !!! Refer to ncdataset and cfdatasetfor additional methods !!!
%
% For more information on the methods use help. For example:
%   >> help ncgeodataset.struct
%   >> doc ncgeodataset
%
% Example 1:
%
%   url='http://geoport.whoi.edu/thredds/dodsC/bathy/gom03_v1_0'
%   nc=ncgeodataset(url)
%   z=nc{'topo'}(500:600,400:500);
%   zg=nc{'topo'}(500:600,400:500).grid;
%   pcolorjw(zg.lon,zg.lat,z);
%
% Example 2:
%
%   url = 'http://dods.mbari.org/cgi-bin/nph-nc/data/ssdsdata/deployments/m1/200810/OS_M1_20081008_TS.nc'
%   geo = ncgeodataset(url)
%   nm = geo.standard_name('sea_water_temperature')
%   v = geo{nm} % works, returns ncgeovariable
%   geo{nm}(1:10,1:2)           % 10x2 {} with matlab indexing
%   geo{nm}(1:10,1:2).data      % 10x2 (works with patch) 
%   geo{nm}(1:10,1:2).grid      % structure with subsetted grid
%   geo.data(nm,[1 1 1 1],[10 2 1 1]) % 10x2 () with netcdf indexing
%   geo.axes(nm)                % Axis names
%   geo.extent(nm)              % geographic extent
%   datestr(geo.timeextent(nm)) % temporal extent
%   tv=geo.gettimevar(nm)       % ncgeovariable of nm's time
%   lonv=geo.getlonvar(nm)      % ncgeovariable of nm's latitude
%   latv=geo.getlatvar(nm)      % ncgeovariable of nm's longitude
%   tnm=geo.gettimename(nm)     % name of nm's time variable
%   lonnm=geo.getlonname(nm)    % name of nm's latitude variable
%   latnm=geo.getlatname(nm)    % name of nm's longitude variable
%   xnm=geo.getxname(nm)        % name of nm's x variable
%   ynm=geo.getyname(nm)        % name of nm's y variable
%
% See also NCTOOLBOX (https://github.com/nctoolbox/nctoolbox
% https://github.com/nctoolbox/nctoolbox), cfdataset, ncdataset, ncgeovariable
% doc ncgeodataset
classdef ncgeodataset < cfdataset
    
    properties (SetAccess = private, GetAccess = private)
        ncvariables
    end
    
    methods
        
        %%
        function obj = ncgeodataset(url)
            % NCGEODATSET  Constructor. Instantiates a NetcdfDataset pointing to the
            % datasource specified by 'url' and uses that as the underlying
            % dataaccess API. When instantiated, the names of all variables
            % are fetched and stored in the 'variables' property. This can be
            % use to open local files, files stored on an HTTP server and
            % OpenDAP URLs.
            obj = obj@cfdataset(url);
            
            % Java hashTable doens't store matlab objects SO we'll use
            % poor man's hash of a n x 2 cell array
            obj.ncvariables = cell(length(obj.variables), 2);
            for i = 1:length(obj.variables)
                obj.ncvariables{i, 1} = obj.variables{i};
            end
        end
        
                %%
        function s = struct(obj, variableName, first, last, stride)
            % NCGEODATASET.STRUCT Retrieve all or a subset of the data for the
            % given variable. The data is returned as a structure containing a
            % variable for the data as well as for each dimension of the
            % data.
            %
            % Usage:
            %   d = cfdataset.struct(variableName)
            %   d = cfdataset.struct(variableName, i1)
            %   d = cfdataset.struct(variableName, i1, i1)
            %   d = cfdataset.struct(variableName, i0, i1, i2)
            %   d = cfdataset.struct(variableName, i0, i1, i2, i3)
            %
            %   If no arguments are provided all the data is returned for the
            %   given variable.
            %
            % Arguments:
            %   variableName = The name of the variable of interest
            %   i1 = index of data for the variables first dimension that
            %   you want returned.
            %   i2 = index of data for the variables second dimension
            %   i3 = index of data for the variables second dimension
            %   i4 = index of data for the variables second dimension
            %   NOTE! first, last, and stride must be matrices the same size as the
            %       matrix returned by NCDATASET.SIZE or SIZE
            %
            % Returns:
            %   The data is returned as a structure containing the actual data for the variable
            %   of interest as well as each coordinate variable
            %
            %
            %   ds = cfdataset('http://dods.mbari.org/cgi-bin/nph-nc/data/ssdsdata/deployments/m1/200810/OS_M1_20081008_TS.nc');
            %   t = ds.struct('TEMP');
            %
            %
            % Note:
            %   The above example is also equivalent to:
            %   ds = cfdataset('http://dods.mbari.org/cgi-bin/nph-nc/data/ssdsdata/deployments/m1/200810/OS_M1_20081008_TS.nc');
            %   v = ds.variable('TEMP')
            %   t = v.data
            %   The reason for providing the 'struct' alternative syntax is
            %   that it may be more familiar to some users.
            v = obj.variable(variableName);
            sz = size(v);
            if (nargin == 2)
              switch length(sz)
                case 1
                  s = v.grid(:);
                  s.(v.name) = v.data(:);
                case 2
                  s = v.grid(:,:);
                  s.(v.name) = v.data(:,:);
                case 3
                  s = v.grid(:,:,:);
                  s.(v.name) = v.data(:,:,:);
                case 4
                  s = v.grid(:,:,:,:);
                  s.(v.name) = v.data(:,:,:,:);
              end
            else
                
                
                % Fill in missing arguments
                % default stride is 1
                if (nargin < 5)
                    stride = ones(1, length(sz));
                end
                
                % Default last is the end
                if (nargin < 4)
                    last = sz;
                end
                
                switch length(sz)
                  case 1
                    s = v.grid(:);
                    s.(v.name) = v.data(first:stride:last);
                  case 2
                    s = v.grid(:,:);
                    s.(v.name) = v.data(first(1):stride(1):last(1),first(2):stride(2):last(2));
                  case 3
                    s = v.grid(:,:,:);
                    s.(v.name) = v.data(first(1):stride(1):last(1),first(2):stride(2):last(2),first(3):stride(3):last(3));
                  case 4
                    s = v.grid(:,:,:,:);
                    s.(v.name) = v.data(first(1):stride(1):last(1),first(2):stride(2):last(2),first(3):stride(3):last(3),...
                      first(4):stride(4):last(4));
                end
                
            end
        end
        
        %%
        function s = grid(obj, variableName, first, last, stride)
            % CFDATASET.GRID Retrieve all or a subset of the coordinate
            % data for the variable. The data is returned as a structure
            % containing a variable for each dimension of the data. The data
            % for the variable is NOT returned ONLY the coordinate data for the
            % variable is returned!!!
            %
            % Usage:
            %   d = cfdataset.grid(variableName)
            %   d = cfdataset.grid(variableName, first)
            %   d = cfdataset.grid(variableName, first, last)
            %   d = cfdataset.grid(variableName, first, last, stride)
            %
            %   If no arguments are provided all the data for the coordinate
            %   variablesis are returned for the given variable.
            %
            % Arguments:
            %   variableName = The name of the variable of interest
            %   first = The first point you want to retrieve (first point idx = 1)
            %   last  = The last point you want to retrive (default is the end of
            %       the data array)
            %   stride = The stride spacing (default is 1)
            %   NOTE! first, last, and stride must be matrices the same size as the
            %       matrix returned by NCDATASET.SIZE or SIZE
            %
            % Returns:
            %   The data is returned as a structure containing the data for all
            %   coordinate variables for the given variable. The data for the
            %   variable itself is NOT returned (see the struct method). THis
            %   is useful if you want to subset based on the coordinate
            %   variables BEFORE fetching over the real data.
            %
            % Example:
            %
            %   ds = cfdataset('http://dods.mbari.org/cgi-bin/nph-nc/data/ssdsdata/deployments/m1/200810/OS_M1_20081008_TS.nc');
            %   td = ds.data('TEMP'); 
            %   tg = ds.grid('TEMP');
            %   ds.attributes('TIME')
            %   plot(datenum(1950, 1, 1) + tg.TIME ,td); datetick('x', 29)
            %
            % Note:
            %   The above example is also equivalent to:
            %   ds = cfdataset('http://dods.mbari.org/cgi-bin/nph-nc/data/ssdsdata/deployments/m1/200810/OS_M1_20081008_TS.nc');
            %   v = ds.variable('TEMP')
            %   tg = v.grid
            %   td = v.data   
            %   The reason for providing the 'grid' alternative syntax is
            %   that it may be more familiar to some users.
            v = obj.variable(variableName);
            sz = size(v);
            
            if (nargin == 2)
              switch length(sz)
                case 1
                 s = v.grid(:);
                case 2
                  s = v.grid(:,:);
                case 3
                  s = v.grid(:,:,:);
                case 4
                  s = v.grid(:,:,:,:);
              end
            else
                
                
                % Fill in missing arguments
                % default stride is 1
                if (nargin < 5)
                    stride = ones(1, length(sz));
                end
                
                % Default last is the end
                if (nargin < 4)
                    last = sz;
                end
                
                switch length(sz)
                  case 1
                    s = v.grid(first:stride:last);
                  case 2
                    s = v.grid(first(1):stride(1):last(1),first(2):stride(2):last(2));
                  case 3
                    s = v.grid(first(1):stride(1):last(1),first(2):stride(2):last(2),first(3):stride(3):last(3));
                  case 4
                    s = v.grid(first(1):stride(1):last(1),first(2):stride(2):last(2),first(3):stride(3):last(3),...
                      first(4):stride(4):last(4));
                end
                
            end
        end
        
        
        %%
        function ax = axes(obj, variableName)
            % NCGEODATASET.AXES(varname)  fetches the names of the axes of varname. 
            % ---- Attempt to fetch the variables representing the axes
            % for the variable of interest. We'll try the CF
            % conventions first and if that's not available we'll try
            % COARDS.
            % See also ncdataset.axes ncvariable.axes ncgeovariable.axes
            
            att = obj.attributes(variableName);
            coordinates = value4key(att, 'coordinates');
            
            if ~isempty(coordinates)
                % ---- Look for CF 'coordinates' attribute
                
                % Parse the string into white space delimited parts
                jls = java.lang.String(coordinates);
                p = jls.split(' ');                   % java.lang.String[]
                axesVariableNames = cell(size(p));    % cell version of p
                for i = 1:length(p)
                    axesVariableNames{i} = char(p(i));
                end
                
            else
                % ---- Look for COARDS conventions. If any coordinate
                % dimensions are missing we don't bother looking any
                % up.
                ncd = ncdataset(obj);
                axesVariableNames = ncd.axes(variableName);
                if ~isempty(axesVariableNames)
                    for i = 1:length(axesVariableNames)
                        if isempty(axesVariableNames{i})
                            axesVariableNames = {};
                            break;
                        end
                    end
                end
                
            end
            
            % Testing combining the axes from both the variable axes and dataset axes
            % maybe temporary stop gap...?
            for i = 1:length(axesVariableNames)
                axesVariables{i,1} = axesVariableNames{i};
            end
            
            if ~exist('axes', 'var')
                %                 try
                ncd = ncdataset(obj);
                dsaxes = ncd.axes(variableName);
                alreadythere = 0;
                for i = 1:length(dsaxes)
                    if ~isempty(dsaxes{i})
                        for j = 1:length(axesVariables)
                            if strcmp(axesVariables{j}, dsaxes{i})
                                alreadythere = 1;
                            end
                        end
                        if ~alreadythere
                            axesVariables{length(axesVariables)+1,1} = dsaxes{i};
                        end
                    end
                end
                
                
                %                     v = ncgeovariable(obj, variableName, axesVariables);
                if ~exist('axesVariables', 'var')
                    ax = {};
                else
                    ax = axesVariables;
                end
                
                
                %                     if ~isempty(v)
                %                         for i = 1:length(obj.variables)
                %                             if strcmp(obj.ncvariables{i, 1}, variableName)
                %                                 obj.ncvariables{i, 2} = v;
                %                                 break;
                %                             end
                %                         end
                %                     end
                %                 catch
                %                     warning('NCGEODATASET:GEOVARIABLE', 'The netcdf-java cdm contains no coordinate information associated with the variable. Returning ncvariable instead of ncgeovariable object. (Methods that rely on coordinate information like ''grid'' or ''geosubset'' are not available.');
                % %                     v = ncvariable(obj, variableName);
                %                 end
            else
                
                %                 v = ncgeovariable(obj, variableName, axes);
            end
        end
       
        
        function v = geovariable(obj, variableName, axes)
            % NCGEODATASET.GEOVARIABLE Returns an ncgeovariable object that provides
            % advanced access to the data contained within that variable based on geo-
            % graphically located data.
            %
            % Usage:
            %    v = ncgeodataset.geovariable(variableName)
            %
            % Arguments:
            %    variableName = A string name of the variable you want to
            %    retrieve. you can use cfdataset.variables to obtain a list
            %    of all variables available.
            %
            % Returns:
            %
            %    v = an instance of ncgeovariable
            %
            % See also ncgeovariable, ncvariable
            %
            
            % Check to see if we've aready fetched the variable of interest
            v = value4key(obj.ncvariables, variableName);
            if isempty(v)
                
                % ---- Attempt to fetch the variables representing the axes
                % for the variable of interest. We'll try the CF
                % conventions first and if that's not available we'll try
                % COARDS.
                
                att = obj.attributes(variableName);
                coordinates = value4key(att, 'coordinates');
                
                if ~isempty(coordinates)
                    % ---- Look for CF 'coordinates' attribute
                    
                    % Parse the string into white space delimited parts
                    jls = java.lang.String(coordinates);
                    p = jls.split(' ');                   % java.lang.String[]
                    axesVariableNames = cell(size(p));    % cell version of p
                    for i = 1:length(p)
                        axesVariableNames{i} = char(p(i));
                    end
                    
                else
                    % ---- Look for COARDS conventions. If any coordinate
                    % dimensions are missing we don't bother looking any
                    % up.
                    ncd = ncdataset(obj);
                    axesVariableNames = ncd.axes(variableName);
                    if ~isempty(axesVariableNames)
                        for i = 1:length(axesVariableNames)
                            if isempty(axesVariableNames{i})
                                axesVariableNames = {};
                                break;
                            end
                        end
                    end
                    
                end
                
                % Testing combining the axes from both the variable axes and dataset axes
                % maybe temporary stop gap...?
                for i = 1:length(axesVariableNames)
                    axesVariables{i,1} = axesVariableNames{i};
                end
                
                if ~exist('axes', 'var')
                    try
                        ncd = ncdataset(obj);
                        dsaxes = ncd.axes(variableName);
                        alreadythere = 0;
                        for i = 1:length(dsaxes)
                            if ~isempty(dsaxes{i})
                                for j = 1:length(axesVariables)
                                    if strcmp(axesVariables{j}, dsaxes{i})
                                        alreadythere = 1;
                                    end
                                end
                                if ~alreadythere
                                    axesVariables{length(axesVariables)+1,1} = dsaxes{i};
                                end
                            end
                        end
                        
                        
                        v = ncgeovariable(obj, variableName, axesVariables);
                        
                        
                        if ~isempty(v)
                            for i = 1:length(obj.variables)
                                if strcmp(obj.ncvariables{i, 1}, variableName)
                                    obj.ncvariables{i, 2} = v;
                                    break;
                                end
                            end
                        end
                    catch
                        warning('NCGEODATASET:GEOVARIABLE', 'The netcdf-java cdm contains no coordinate information associated with the variable. Returning ncvariable instead of ncgeovariable object. (Methods that rely on coordinate information like ''grid'' or ''geosubset'' are not available.');
                        v = ncvariable(obj, variableName);
                    end
                else
                    
                    v = ncgeovariable(obj, variableName, axes);
                end
                
            end
        end
        
        function e = extent(obj, variableName)
            % NCGEODATASET.extent(variable) - Function to calculate lat/lon bounding box of variable.
            % Usage: ex = nc.extent('temp')
            v = obj.geovariable(variableName);
            e = v.extent;
            
        end
        
        function e = top_bot(obj, variableName)
            % NCGEODATASET.top_bot - Function to calculate lat/lon bounding box of variable.
            % Usage: ex = nc.extent('temp')
            v = obj.geovariable(variableName);
            e = v.top_bot;
            
        end
        
        function te = timeextent(obj, variableName)
            % NCGEODATASET.timeextent - Function to calculate the start and stop times of for a variable.
            % Usage: t = nc.timeextent('salt')
            v = obj.geovariable(variableName);
            te = v.timeextent;
        end
        
        function tv = gettimevar(obj, variableName)
            % NCGEODATASET.gettimevar
            var = obj.geovariable(variableName);
            tv = var.gettimevar();
        end
        
        function tv = getlonvar(obj, variableName)
            % NCGEODATASET.gettimevar
            var = obj.geovariable(variableName);
            tv = var.getlonvar();
        end
        
        function tv = getlatvar(obj, variableName)
            % NCGEODATASET.gettimevar
            var = obj.geovariable(variableName);
            tv = var.getlatvar();
        end
        
        function tv = gettimename(obj, variableName)
            % NCGEODATASET.gettimename
            var = obj.geovariable(variableName);
            tv = var.gettimename();
        end
        
        function tv = getlatname(obj, variableName)
            % NCGEODATASET.gettimename
            var = obj.geovariable(variableName);
            tv = var.getlatname();
        end
        
        function tv = getyname(obj, variableName)
            var = obj.geovariable(variableName);
            tv = var.getyname();
        end
        
        function tv = getlonname(obj, variableName)
            % NCGEODATASET.gettimename
            var = obj.geovariable(variableName);
            tv = var.getlonname();
        end
        
        function tv = getxname(obj, variableName)
            var = obj.geovariable(variableName);
            tv = var.getxname();
        end
        
        function tn = gettimedata(obj, variableName, start, last, stride)
            % NCGEOVARIABLE.gettimedata()
            var = obj.getimevar(variableName);
            tn = var.gettimedata(start, last, stride);
        end
        
        function s = getlondata(obj, variableName, start, last, stride)
            % NCGEOVARIABLE.getlondata()
            var = obj.geovariable(variableName);
            sz = var.size();
            switch length(sz)
                case 1
                    s = var.getlondata(start:stride:last);
                case 2
                    s = var.getlondata(start(1):stride(1):last(1),start(2):stride(2):last(2));
                case 3
                    s = var.getlondata(start(1):stride(1):last(1),start(2):stride(2):last(2),start(3):stride(3):last(3));
                case 4
                    s = var.getlondata(start(1):stride(1):last(1),start(2):stride(2):last(2),start(3):stride(3):last(3),...
                        start(4):stride(4):last(4));
            end
            
        end
        
        function s = getlatdata(obj, variableName, start, last, stride)
            % NCGEOVARIABLE.getlatdata()
            var = obj.geovariable(variableName);
            sz = var.size();
            switch length(sz)
                case 1
                    s = var.getlatdata(start:stride:last);
                case 2
                    s = var.getlatdata(start(1):stride(1):last(1),start(2):stride(2):last(2));
                case 3
                    s = var.getlatdata(start(1):stride(1):last(1),start(2):stride(2):last(2),start(3):stride(3):last(3));
                case 4
                    s = var.getlatdata(start(1):stride(1):last(1),start(2):stride(2):last(2),start(3):stride(3):last(3),...
                        start(4):stride(4):last(4));
            end
        end
        
        
        %% Should be called as ncpoint(nc), etc.
        %
        %         function p = point(obj)
        %             p = ncpoint(obj);
        %         end
        %
        %         function r = rgrid(obj)
        %             r = ncrgrid(obj);
        %         end
        %
        %         function c = cgrid(obj)
        %             c = nccgrid(obj);
        %         end
        %
        %         function sg = sgrid(obj)
        %             sg = ncsgrid(obj);
        %         end
        %
        %         function ug = ugrid(obj)
        %             ug = ncugrid(obj);
        %         end
        %%
        % CFDATASET.NUMEL Overridden function required for supporting
        % SUBSREF
        function result = numel(varargin)
            % cfdataset/numel -- Overloaded NUMEL.
            result = 1;
        end
        
        function B = subsref(obj,s)
            if isa(s(1).subs, 'cell')
                s(1).subs = s(1).subs{1};
            end
            
            switch s(1).type
                case '.'
                    if length(s) < 2
                        B = builtin('subsref',obj,s);
                    elseif length(s) == 2
                        if s(2).type == '.'
                            a = substruct('.',s(1).subs);
                            A = builtin('subsref',obj,a);
                            %                         b = substruct('.',s(2).subs,'()',s(2).subs);
                            B = builtin('subsref',A,s(2:end));
                        else
                            g = substruct('.',s(1).subs,'()',s(2).subs);
                            % g.type = '()';
                            % g.subs = s(i).subs;
                            
                            B = builtin('subsref',obj,g);
                        end
                    elseif length(s)==3
                        switch s(2).type
                            case '.'
                                a = substruct('.',s(1).subs);
                                A = builtin('subsref',obj,a);
                                B = builtin('subsref',A,s(2:end));
                            case '()'
                                B = builtin('subsref',obj,s);
                                
                        end
                    else
                        error(['NCTOOLBOX:ncgeodataset:subsref'], ...
                            'Unexpected reference notation or method call')
                    end
                    
                case '()'
                    error('NCTOOLBOX:ncgeodataset:subsref',...
                        'Call with "()" as first type unsupported at this time')
                    
                case '{}'
                    warning('off', 'all');
                    echo off
                    v = obj.geovariable(s(1).subs);
                    echo on
                    if length(s) == 1
                        B = v;
                    elseif length(s) == 2
                        echo off
                        B = v.data(s(2).subs{:});
                        echo on
                    elseif length(s) == 3
                        echo off
                        switch s(3).subs
                            case 'data'
                                B = squeeze(v.data(s(2).subs{:})); % Adding squeeze of results related to Issue 23
                            case 'grid'
                                A = v.grid_interop(s(2).subs{:});
                                try %Filtering added for njtbx similar results, entire syntax will be deprecated in the future
                                    B.time = squeeze(A.time); % Adding squeeze of results related to Issue 23
                                catch me
                                end
                                try
                                    B.lon = squeeze(A.lon); % Adding squeeze of results related to Issue 23
                                catch me
                                end
                                try
                                    B.lat = squeeze(A.lat); % Adding squeeze of results related to Issue 23
                                catch me
                                end
                                try
                                    B.z = squeeze(A.z); % Adding squeeze of results related to Issue 23
                                catch me
                                end
                                try
                                    B.x = squeeze(A.x); % Adding squeeze of results related to Issue 23
                                catch me
                                end
                                try
                                    B.y = squeeze(A.y); % Adding squeeze of results related to Issue 23
                                catch me
                                end
                        end
                        echo on
                    else
                        B = v;
                        warning('on', 'all');
                    end
                    
            end
        end
        
    end  % methods end
    
    
end % class end
