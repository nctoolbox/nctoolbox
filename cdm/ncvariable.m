% NCVARIABLE Provide advanced access to variables and their related
% dimensions.
%
% NCVARIABLE is used to retrieve data for a given variable as well as the 
% variables associated coordinate dimensions. Normally, you would retrive
% it using CFDATASET.VARIABLE 
%
% See also CFDATASET
classdef ncvariable < handle
    
    properties (SetAccess = private)
        dataset          % ncdataset instance
    end
    
    properties (Dependent = true)
        name            % The string variable name that this object represents
        axes
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
                ex = MException('NCVARIABLE:ncvariable', 'Invalid dataset was specified');
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
            v = char(obj.variable.getName());
        end
        
        %%
        function n = size(obj)
            % NCVARIABLE.SIZE returns the size of the variable, including
            % it's singleton dimensions
            n = obj.dataset.size(obj.name);
        end
        
        %%
        function d = data(obj, first, last, stride)
            % NCVARIABLE.DATA Retrieve all or a subset of the data for the
            % variable. The data is returned as a nested structure with the
            % 2 substructures of:
            %   metadata
            %       variable = the name of the variable
            %       coordinates = A cell array of the names of the coordinate variables in the correct order
            %   data = the actual data for the variable and each coordinate variable
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
            %   The data is returned as a structure containing the actual data for the variable 
            %   of interest as well as each coordinate variable
            %
            % Example:
            %
            %   ds = cfdataset('http://dods.mbari.org/cgi-bin/nph-nc/data/ssdsdata/deployments/m1/200810/OS_M1_20081008_TS.nc');
            %   v = ds.variable('TEMP');
            %   t = v.data([1 1 1 1], [10 2 1 1]);
            %
            
            if (nargin == 1)
                d = alldata(obj);
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
                
                d = somedata(obj, first, last, stride);
            end
        end
        
        %% see http://www.mathworks.com/access/helpdesk/help/techdoc/index.html?/access/helpdesk/help/techdoc/ref/subsref.html&http://www.mathworks.com/access/helpdesk/help/techdoc/matlab_oop/br02znk-1.html
        function sref = subsref(obj, s)
            % obj(i) is equivalent to obj.Data(i)
            switch s(1).type
                % Use the built-in subsref for dot notation
                case '.'
                    sref = builtin('subsref',obj,s);
                case '()'
                    if length(s)<2
                        % Note that obj.Data is passed to subsref
                        sref = builtin('subsref',obj.Data,s);
                        return
                    else
                        sref = builtin('subsref',obj,s);
                    end
                    % No support for indexing using '{}'
                case '{}'
                    error('NCVARIABLE:subsref',...
                        'Not a supported subscripted reference')
            end
        end
    end
    
    methods (Access = protected)
        %%
        function data = alldata(obj)
            
            % ---- Step 2: Add the data
            name = char(obj.variable.getName());
            data.(name) = obj.dataset.data(name);
            
            for i = 1:length(obj.axesVariables)
                name = char(obj.axesVariables{i}.getName());
                data.(name) = obj.dataset.data(name);
            end
        end
        
        %%
        function data = somedata(obj, first, last, stride)
            
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
            name = char(obj.variable.getName());
            data.(name) = obj.dataset.data(name, first, last, stride);

            % ---- Step 3: Add the data for each axes variable
            for i = 1:length(obj.axesVariables)
                name = char(obj.axesVariables{i}.getName());
                
                % ---- Step 4: figure out how to subset the data properly
                vs = obj.dataset.size(name);
                if (length(vs) == length(s)) 
                    %% case: sizes are the same
                    if isequal(vs, s)
                        vFirst = first;
                        vLast = last;
                        vStride = stride;
                    else
                        me = MException('NCVARIABLE:somedata', ['The data size of the coordinate variable,' ...
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
                        me = MException('NCVARIABLE:somedata', ['The data size of the coordinate variable,' ...
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
                                me = MException('NCVARIABLE:somedata', ['The data size of the coordinate variable,' ...
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
            end
        end
    end
    
end