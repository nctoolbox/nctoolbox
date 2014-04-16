% NCVARIABLE Provide advanced access to variables and their related
% dimensions.
%
% NCVARIABLE is used to retrieve data for a given variable as well as the
% variables associated coordinate dimensions. Normally, you would retrive
% it using CFDATASET.VARIABLE
%
% Example of use:
%  ds = cfdataset('http://dods.mbari.org/cgi-bin/nph-nc/data/ssdsdata/deployments/m1/200810/OS_M1_20081008_TS.nc');
%  v = ds.variable('TEMP');
%
%  % Look at properties
%  v.name
%  v.axes
%
%  % Data access example #1
%  temp = v.data([1 1 1 1], [100 5 1 1]);
%
%  % Data access example #2
%  t_end = v.size;
%  t_start = t_end ./ t_end;
%  t_stride = t_start;
%  t_stride(1)=10;
%  t = v.data(t_start, t_end, t_stride);
%
%
% See also CFDATASET, SIZE, DATA
% NCTOOLBOX (https://github.com/nctoolbox/nctoolbox)
classdef ncvariable < handle
    
    properties (SetAccess = private)
        dataset          % ncdataset instance
    end
    
    properties (Dependent = true)
        name            % The string variable name that this object represents
        axes2           %  ????
        axes            % the coordinate variables associated with the object
        attributes      % The attributes associated with the object.
    end
    
    properties (SetAccess = private, GetAccess = protected)
        variable        % ucar.nc2.Variable instance. Represents the data
        axesVariables    % ucar.nc2.Variable instance. Represents the data.
    end
    
    methods
        
        %%
        function obj = ncvariable(src, variableName, axesVariableNames)
            % NCVARIABLE.NCVARIABLE  Constructor.
            %
            % Use as:
            %    v = ncvariable(src, variableName)
            %    v = ncvariable(src, variableName, axesVariableNames)
            %
            if ischar(src)
                obj.dataset = ncdataset(src);  % src is a string URL/File
            elseif isa(src, 'ncdataset')
                obj.dataset = src;             % src is a ncdataset
            else
                ex = MException('NCTOOLBOX:ncvariable', 'Invalid dataset was specified');
                ex.throw;
            end
            
            obj.variable = obj.dataset.netcdf.findVariable(variableName);
            
            if nargin == 3
                obj.axesVariables = cell(size(axesVariableNames));
                for i = 1:length(axesVariableNames)
                    obj.axesVariables{i} = obj.dataset.netcdf.findVariable(axesVariableNames{i});
                end
            else
                obj.axesVariables = {};
            end
            
            
        end
        
        function a = get.attributes(obj)
            % NCVARIABLE.ATTRIBUTES returns the attributes for the variable.
            a = obj.dataset.attributes(obj.name);
        end
        
        %%
        function a = get.axes(obj)
            % NCVARIABLE.AXES Returns the names of the coordinate axes
            a = cell(size(obj.axesVariables));
            for i = 1:length(obj.axesVariables)
                name = char(obj.axesVariables{i}.getName());
                a{i} = name;
            end
        end
        
        %%
        function v = get.name(obj)
            % NCVARIABLE.NAME Provides dynamic access to the underlying
            % netcdf datasets variable name
            %
            % Example:
            %   vname=v.name % returns the name of the variable
            %
            v = char(obj.variable.getName());
        end
        
        %%
        function n = size(obj)
            % NCVARIABLE.SIZE returns the size of the variable, including
            % its singleton dimensions
            n = obj.dataset.size(obj.name);
        end
        
        function val = attribute(obj, key)
            % NCVARIABLE.ATTRIBUTE returns the value a global attribute specified by its key or the
            % variable attribute specified by key and variable.
            %
            % Use as:
            %   a = ncvariable.attribute('title')
            %   a = ncvariable.attribute(key)
            %
            %
            % Inputs:
            %   key = The name of the attribute field like 'title' or 'units'...
            %
            % Return:
            %   The value associated with the attribute field corresponding to key.
            val = obj.dataset.attribute(obj.name, key);
        end
        
        %%
        function d = data(obj, first, last, stride)
            % NCVARIABLE.DATA Retrieve all or a subset of the data for the
            % variable. The data is returned as a structure containing a
            % variable for the data as well as for each dimension of the
            % data.
            %
            % Usage:
            %   d = ncvariable.data
            %   d = ncvariable.data(first)
            %   d = ncvariable.data(first, last)
            %   d = ncvariable.data(first, last, stride)
            %
            %   If no arguments are provided all the data is returned.
            %
            % Arguments:
            %   first = The first point you want to retrieve (first point idx = 1)
            %   last  = The last point you want to retrive (default is the end of
            %       the data array)
            %   stride = The stride spacing (default is 1)
            %   NOTE! first, last, and stride must be matrices the same size as the
            %       matrix returned by NCDATASET.SIZE or SIZE
            %
            % Returns:
            %   The data is returned as a variable containing the actual data for the netcdf variable
            %
            % Example:
            %
            %   ds = cfdataset('http://dods.mbari.org/cgi-bin/nph-nc/data/ssdsdata/deployments/m1/200810/OS_M1_20081008_TS.nc');
            %   v = ds.variable('TEMP');
            %   t = v.data([1 1 1 1], [10 2 1 1]);
            %
            
            if (nargin == 1)
                if isempty(obj.size) % check if variable obj is dimensionless (not cf compliant.....)
                    try
                        d = alldata(obj,1);
                    catch me
                        ex = MException('NCTOOLBOX:ncvariable:data', ['Failed to open ' url]);
                        ex = ex.addCause(me);
                        ex.throw;
                    end
                else
                    d = alldata(obj, 1);
                end
            else
                s = obj.size;
                
                % Fill in missing arguments
                % default stride is 1
                if (nargin < 4)
                    stride = ones(1, length(s));
                end
                
                % Default last is the end
                if (nargin < 3)
                    last = s;
                end
                
                d = somedata(obj, 1, first, last, stride);
            end
        end
        
        function g = grid(obj, first, last, stride)
            % NCVARIABLE.GRID Retrieve all or a subset of the coordinate
            % data for the variable. The data is returned as a structure
            % containing a variable for each dimension of the data.
            %
            % Usage:
            %   d = ncvariable.grid
            %   d = ncvariable.grid(first)
            %   d = ncvariable.grid(first, last)
            %   d = ncvariable.grid(first, last, stride)
            %
            %   If no arguments are provided all the data is returned.
            %
            % Arguments:
            %   first = The first point you want to retrieve (first point idx = 1)
            %   last  = The last point you want to retrive (default is the end of
            %       the data array)
            %   stride = The stride spacing (default is 1)
            %   NOTE! first, last, and stride must be matrices the same size as the
            %       matrix returned by NCDATASET.SIZE or SIZE
            %
            % Returns:
            %   The data is returned as a structure containing the actual data for the variable
            %   of interest as well as each coordinate variable
            %
            % Example:
            %
            %   ds = cfdataset('http://dods.mbari.org/cgi-bin/nph-nc/data/ssdsdata/deployments/m1/200810/OS_M1_20081008_TS.nc');
            %   v = ds.variable('TEMP');
            %   v.size  % 9043x11x1x1
            %   g = v.grid
            %   t = v.data([1 1 1 1], [10 2 1 1]); % is not like
            %
            
            if (nargin == 1)
                g = alldata(obj, 0);
            else % all this stuff needs logic work to function as expected
                s = obj.size;
                
                % Fill in missing arguments
                % default stride is 1
                if (nargin < 4)
                    stride = ones(1, length(s));
                end
                
                % Default last is the end
                if (nargin < 3)
                    last = s;
                end
                
                g = somedata(obj, 0, first, last, stride);
            end
        end
       
        %%
        function e = end(obj, k, n)
            % NCVARIABLE.END the last index in an indexing
            % expression.
            % v.end(2) % returns the last index of the second dimension.
            % e.g.: elevation.data(end-3:end,1,1)
            n = obj.dataset.size(obj.name);
            e = n(k);
        end % Added to deal with end indexing functionality,
        % otherwise the indexing arugment is ignored.
        
        %%
        function sref = subsref(obj,s)
            %                disp(s(2).subs{3})
            % SUBSREF parses an object name for .name or ()
            switch s(1).type
                % Use the built-in subsref for dot notation
                case '.'
                    switch s(1).subs
                        case 'data'
                            %% .data
                            nums = size(obj);
                            if ~isempty(nums)
                                switch length(s)
                                    case 1
                                        sref = obj;
                                    case 2
                                        idx = s(2).subs;
                                        nidx = length(idx);
                                        % Fill in missing arguments
                                        % default stride is 1
                                        if (nidx < 3)
                                            stride = ones(1, length(nums));
                                        else
                                            stride = idx{3};
                                        end

                                        % Default last is the end
                                        if (nidx < 2)
                                            last = nums;
                                        else
                                            last = idx{2};
                                        end
                                        
                                        first = idx{1};
                                        sref = obj.somedata(1, first, last, stride);
                                end
                                
                            else
                                sref = obj.data;
                                warning('NCTOOLBOX:ncvariable:subsref', ...
                                    ['Variable "' obj.name '" has no netcdf dimension associated with it. Errors may result from non CF compliant files.'])
                            end
                        
                        case 'grid'
                            %% .grid
                            nums = size(obj);
                            if ~isempty(nums)
                                switch length(s)
                                    case 1
                                        sref = obj;
                                    case 2
                                        idx = s(2).subs;
                                        nidx = length(idx);
                                        % Fill in missing arguments
                                        % default stride is 1
                                        if (nidx < 3)
                                            stride = ones(1, length(nums));
                                        else
                                            stride = idx{3};
                                        end

                                        % Default last is the end
                                        if (nidx < 2)
                                            last = nums;
                                        else
                                            last = idx{2};
                                        end
                                        first = idx{1};
                                        sref = obj.grid(first, last, stride);
                                end
                                
                            else
                                warning('NCTOOLBOX:ncvariable:subsref', ...
                                    ['Variable "' name '" has no netcdf dimension associated with it. Errors may result from non CF compliant files.'])
                            end
                        otherwise
                            sref = builtin('subsref',obj,s);
                    end
                case '()'
                    if length(s)<2
                        % Note that obj.Data is passed to subsref
                        sref = builtin('subsref',obj.data,s);
                    else
                        sref = builtin('subsref',obj,s);
                    end
                    % No support for indexing using '{}'
                case '{}'
                    error('NCTOOLBOX:ncvariable:subsref', ...
                        'Not a supported subscripted reference, "{}" are not permitted to call variable object methods');
            end
        end
        
        
    end
    
    methods (Access = protected)
        
        
        %%
        function data = alldata(obj, withData)
            % NCVARIABLE.ALLDATA -- extract the data or the axes
            % v.alldata(1) % returns the data
            % v.alldata()  % returns a structure with the axesVariables
            %              % and their data
            if withData == 1
                name = char(obj.variable.getName());
                data = obj.dataset.data(name);
            end
            
            if withData == 0
                for i = 1:length(obj.axesVariables)
                    name = char(obj.axesVariables{i}.getName());
                    data.(name) = obj.dataset.data(name);
                end
            end
        end
        
        %%
        function data = somedata(obj, withData, first, last, stride)
            % NCGEOVARIABLE.SOMEDATA -- extract a hyperslab of data
            % SOMEDATA can extract a subset of data from an NCVARIABLE or
            % the subset of the corresponding axis variables.
            % 
            % Example:
            %  first=ones(length(v.size));
            %  last=v.size;
            %       stride=first;
            % vd = v.somedata(1,first,last,stride); % returns data
            % vg = v.somedata(0,first,last,stride); % returns struct w/grid
            % 
            % See also NCVARIABLE, NCVARIABLE.ALLDATA,
            % NCVARIABLE.DATA, NCVARIABLE.GRID
            s = obj.dataset.size(obj.name);
            
            % Fill in missing arguments
            % default stride is 1
            if (nargin < 5)
                stride = ones(1, length(s));
            end
            
            % Default last is the end
            if (nargin < 4)
                last = s;
            end
            
            % ---- Step 2: Add the data for the variable of interest
            if withData == 1
                name = char(obj.variable.getName());
                data = obj.dataset.data(name, first, last, stride);
            end
            
            % ---- Step 3: Add the data for each axes variable
            if withData == 0
                for i = 1:length(obj.axesVariables)
                    name = char(obj.axesVariables{i}.getName());
                    %                     type = char(obj.axesVariables{i}.getAxisType());
                    % ---- Step 4: figure out how to subset the data properly
                    vs = obj.dataset.size(name);
                    if numel(vs) > 0 % Added to solve work around somedata calls that involve variables with
                        % no netcdf dim. (This will be frequent in some OOI files.)
                        if (length(vs) == length(s))
                            %% case: sizes are the same
                            if isequal(vs, s)
                                vFirst = first;
                                vLast = last;
                                vStride = stride;
                            else
                                me = MException('NCTOOLBOX:ncvariable:somedata', ...
                                    ['The data size of the coordinate variable,' ...
                                    name ', does not fit the size of ' obj.name]);
                                me.throw;
                            end
                            
                        elseif length(vs) == 1
                            %% case: singleton dimension. Find side of data with
                            % the same length
                            
                            % TODO: the following line  will give bogus results if
                            % the data has 2 dimensions of the same length
                            dim = find(s == vs, 1);
                            if ~isempty(dim)
                                vFirst = first(dim);
                                vLast = last(dim);
                                vStride = stride(dim);
                            else
                                me = MException('NCTOOLBOX:ncvariable:somedata', ...
                                    ['The data size of the coordinate variable,' ...
                                    name ', does not fit the size of ' obj.name]);
                                me.throw;
                            end
                            
                        else
                            %% case: variable is coordiantes. Look for size
                            % TODO this is a lame implementation.
                            dim = find(s == vs(1), 1);
                            if ~isempty(dim)
                                for j = 2:length(vs)
                                    if vs(j) ~= s(dim + j - 1)
                                        me = MException('NCTOOLBOX:ncvariable:somedata', ...
                                            ['The data size of the coordinate variable,' ...
                                            name ', does not fit the size of ' obj.name]);
                                        me.throw;
                                    end
                                end
                                k = dim:dim + length(vs) - 1;
                                vFirst = first(k);
                                vLast = last(k);
                                vStride = stride(k);
                            end
                            
                        end
                        data.(name) = obj.dataset.data(name, vFirst, vLast, vStride);
                    else
                        data.(name) = obj.dataset.data(name);
                    end
                end
            end
        end
    end
    
end
