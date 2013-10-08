% CFFEATURE
%
% NCTOOLBOX (http://github.com/nctoolbox/nctoolbox)
classdef cffeature < cfdataset
    
    properties (SetAccess = private, GetAccess = private)
        ncvariables
        points
        feature_collection
        feature_type
        time
        lat
        lon
        z
    end
    
    methods
        
        %%
        function obj = cffeature(url)
            % CFFEATURE  Constructor. Instantiates a NetcdfDataset pointing to the
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
            
            % Keep a reference to the feature wrapped dataset
            wrapped = wrapfeature(obj.netcdf);
            wrapped = wrapped.getPointFeatureCollectionList();
            feature = wrapped.get(0);
            obj.feature_type = feature.getCollectionFeatureType();
            obj.feature_collection = feature;
        end
        
        %%
        function e = extent(obj, variableName)
            % CFFEATURE.extent - Function to calculate lat/lon bounding box of variable.
            % Usage: ex = nc.extent('temp')
            ft = obj.feature_type;
            if strcmp(ft, 'PROFILE') | strcmp(ft, 'TRAJECTORY') | strcmp(ft, 'SECTION') | strcmp(ft, 'POINT')
                if isempty(obj.points)
                    setpoints(obj)
                end
                lat = obj.getlatdata();
                lon = obj.getlondata();
                e.lon = [min(lon), max(lon)];
                e.lat = [min(lat), max(lat)];
            else
                javabbox = obj.feature_collection.getBoundingBox();
                e.lon = [javabbox.getLonMin, javabbox.getLonMax];
                e.lat = [javabbox.getLatMin, javabbox.getLatMax];
            end
        end
        
        function te = timeextent(obj, variableName)
            % CFFEATURE.timeextent - Function to calculate the start and stop times of for a variable.
            % Usage: t = nc.timeextent('salt')
            t = obj.gettimedata();
            te = [min(t), max(t)];
        end
        
        function locs = positions(obj)
            locs.lon = obj.getlondata();
            locs.lat = obj.getlatdata();
        end
        
        function d = data(obj, variableName)
            for i = 1:length(obj.points)
               d(i) = str2num(obj.points(i).getData().getArray(varname).toString());
            end
        end
        
        function g = grid()
            g.lat = obj.getlatdata();
            g.lon = obj.getlondata();
            g.z = obj.getelevationdata();
            g.time = obj.gettimedata();
        end
        
        function s = geosubset()
            % TODO Implement
        end
        
        function tn = gettimedata(obj, variableName)
            % CFFEATURE.gettimedata()
            if isempty(obj.time)
                if isempty(obj.points)
                    setpoints(obj)
                end
%                 for i = 1:length(obj.points)
%                     ds{i} = char(obj.points(i).getNominalTimeAsCalendarDate());
%                 end
%                 obj.time = datenum(ds, 'yyyy-mm-ddTHH:MM:SSZ');
            end
            tn = obj.time;
        end
        
        function ln = getlondata(obj, variableName)
            % CFFEATURE.getlondata()
            if isempty(obj.lat)
                if isempty(obj.points)
                   setpoints(obj); 
                end
%                 obj.lat = [];
%                 for i = 1:length(obj.points)
%                     loc = obj.points(i).getLocation();
%                     obj.lat(i) = loc.getLatitude();
%                     obj.lon(i) = loc.getLongitude();
%                 end
            end
            ln = obj.lon;
        end
        
        function ln = getlatdata(obj, variableName)
            % CFFEATURE.gelatdata()
            if isempty(obj.lat)
                if isempty(obj.points)
                   setpoints(obj); 
                end
%                 obj.lon = [];
%                 for i = 1:length(obj.points)
%                     loc = obj.points(i).getLocation();
%                     obj.lat(i) = loc.getLatitude();
%                     obj.lon(i) = loc.getLongitude();
%                 end
            end
            ln = obj.lat;
        end
        
        function z = getelevationdata()
            if isempty(obj.z)
                if isempty(obj.points)
                    setpoints(obj);
                end
%                 for i = 1:length(obj.points)
%                     loc = obj.points(i).getLocation();
%                     obj.z(i) = loc.getAltitude();
%                 end
            end
            z = obj.z;
        end
        
        function s = setpoints(obj)
            if ~isempty(strfind(obj.feature_type, 'STATION'))
                stations = obj.feature_collection.getStations();
                obj.points = iterstations(stations);
            else
                obj.points = iterfeature(obj.feature_collection);
            end
            for i = 1:length(obj.points)
                ds{i} = char(obj.points(i).getNominalTimeAsCalendarDate());
                loc = obj.points(i).getLocation();
                obj.lat(i) = loc.getLatitude();
                obj.lon(i) = loc.getLongitude();
                obj.z(i) = loc.getAltitude();
            end
            obj.time = datenum(ds, 'yyyy-mm-ddTHH:MM:SSZ');
        end
        
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
                        error(['NCTOOLBOX:CFFEATURE:subsref'], ...
                            'Unexpected reference notation or method call')
                    end
                    
                case '()'
                    error('NCTOOLBOX:CFFEATURE:subsref',...
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
                                B = squeeze(v.data(s(2).subs)); % Adding squeeze of results related to Issue 23
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
