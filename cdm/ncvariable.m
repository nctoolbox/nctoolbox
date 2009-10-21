classdef ncvariable < handle
    
    properties (SetAccess = private)
        dataset          % ncdataset instance
    end
    
    properties (Dependent = true)
        name         % The string variable name that this object represents
    end
    
    properties (SetAccess = private, GetAccess = protected)
        variable     % ucar.nc2.Variable instance. Represents the data
        axesVariables    % ucar.nc2.Variable instance. Represents the data.
    end
    
    methods
        function obj = ncvariable(src, variableName, axesVariableNames)
            % NCVARIABLE  Constructor.
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
                ex = MException('MBARI:NCVARIABLE', 'Invalid dataset was specified');
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
        
        function v = get.name(obj)
            v = char(obj.variable.getName());
        end
        
        function n = size(obj)
            n = obj.dataset.size(obj.name);
        end
        
        function d = data(obj, first, last, stride)
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
                
                d = alldata(obj, first, last, stride);
            end
        end
        
        % see http://www.mathworks.com/access/helpdesk/help/techdoc/index.html?/access/helpdesk/help/techdoc/ref/subsref.html&http://www.mathworks.com/access/helpdesk/help/techdoc/matlab_oop/br02znk-1.html
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
                    error('MYDataClass:subsref',...
                        'Not a supported subscripted reference')
            end
        end
    end
    
    methods (Access = protected)
        function d = alldata(obj)
            % ---- Step 1: Specify the name of the data
            d.metadata.variable = obj.name;
            
            % ---- Step 2: Add the data
            name = char(obj.variable.getName());
            d.data.(name) = obj.dataset.data(name);
            
            % ---- Step 3: Add the coordinate dimensions in order as strings
            % ---- Step 4: Add the cooridinate data
            d.metadata.coordinates = cell(size(obj.axesVariables));
            for i = 1:length(obj.axesVariables)
                name = char(obj.axesVariables{i}.getName());
                d.metadata.coordinates{i} = name;
                d.data.(name) = obj.dataset.data(name);
            end
        end
        
        function d = somedata(obj, first, last, stride)
            
            s = obj.size(obj.name);
            
            % Fill in missing arguments
            % default stride is 1
            if (nargin < 5)
                stride = ones(1, length(s));
            end
            
            % Default last is the end
            if (nargin < 4)
                last = s;
            end
            
            % ---- Step 1: Specify the name of the data
            d.metadata.variable = obj.name;
            
            % ---- Step 2: Add the data
            name = char(obj.variable.getName());
            d.data.(name) = obj.dataset.data(name, first, last, stride);
            
            % ---- Step 3: Add the coordinate dimensions in order as strings
            % ---- Step 4: Add the cooridinate data
            d.metadata.coordinates = cell(size(obj.axesVariables));
            for i = 1:length(obj.axesVariables)
                name = char(obj.axesVariables{i}.getName());
                d.metadata.coordinates{i} = name;
                
                % ---- Step 5: figure out how to subset the data properly
                
                
                d.data.(name) = obj.dataset.data(name);
            end
        end
    end
    
end