% NCGEODATASET
%
% NCTOOLBOX (http://code.google.com/p/nctoolbox)
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
        function ax = axes(obj, variableName)
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
                ax = axesVariables;
                
                
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
            % NCGEODATASET.VARIABLE Returns an ncgeovariable object that provides
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
            % NCGEODATASET.extent - Function to calculate lat/lon bounding box of variable.
            % Usage: ex = nc.extent('temp')
            v = obj.geovariable(variableName);
            e = v.extent;
            
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
        
        function tv = getlonname(obj, variableName)
            % NCGEODATASET.gettimename
            var = obj.geovariable(variableName);
            tv = var.getlonname();
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
        
        function ln = getlatdata(obj, variableName, start, last, stride)
            % NCGEOVARIABLE.gelatdata()
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
