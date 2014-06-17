classdef cdm < handle

    properties(SetAccess = private)
        location
        netcdf
        variables
    end

    properties(SetAccess = protected, GetAccess = ?cdmvariable)
        dataset % A CFDATASET
    end

    methods

        function obj = cdm(url)
            if isa(url, 'ncdataset')
                obj.dataset = cfdataset(url.location);
            elseif isa(url, 'char')
                obj.dataset = cfdataset(url);
            end
            obj.location = obj.dataset.location;
            obj.netcdf = obj.dataset.netcdf;
            obj.variables = obj.dataset.variables;
        end

        % Delegate to contained CFDATASET
        function v = size(obj, variable)
            v = obj.dataset.size(variable);
        end

        function v = axes(obj, variable)
            v = obj.dataset.axes(variable);
        end

        function v = dimensions(obj, variable)
            v = obj.dataset.dimensions(variable);
        end

        function v = attributes(obj, variable)
            v = obj.dataset.attributes(variable);
        end

        function v = attribute(obj, varargin)
            v = obj.dataset.attribute(varargin);
        end

        function v = metadata(obj)
            v = obj.dataset.metadata;
        end

        function save(obj, filename)
            v = obj.save(filename);
        end

        function v = time (obj, variable, data)
            v = obj.time(variable, data);
        end

        function v = ncml(obj)
            v = obj.dataset.ncml;
        end
        
        % Methods to be converted to matlab style indexing from netcdf style indexing
        % e.g. function d = data(obj, variable, first, last, stride)

        function v = variable(obj, variableName)
            v = cdmvariable(obj, variableName); 
        end

        function v = data(obj, variableName, varargin)
            [first, last, stride] = obj.toncindex(variableName, varargin);
            v = obj.dataset.data(variableName, first, last, stride);
        end

        function v = struct(obj, variableName, varargin)
            [first, last, stride] = obj.toncindex(variableName, varargin);
            v = obj.dataset.struct(variableName, first, last, stride);
        end

        function v = grid(obj, variableName, varargin)
            [first, last, stride] = obj.toncindex(variableName, varargin);
            v = obj.dataset.grid(variableName, first, last, stride);
        end

        function v = standard_name(obj, standardName)
            v = obj.dataset.standard_name(standardName);
        end
        

    end

    methods (Access = protected)
        function [first, last, stride] = toncindex(obj, variableName, a)
            osize = double(size(obj.variable(variableName)));
            [~, vc] = size(a);
            if vc == 0   % handle .mdata and mdata() cases
                a = {{':'}};
            end
            [first, last, stride] = indexing(a, osize);
        end
        
    end

end