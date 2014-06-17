% CFVARIABLE Provide advanced access to variables and their related
% dimensions.
%
% CFVARIABLE is used to retrieve data for a given variable as well as the
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
classdef cdmvariable < handle

    properties (SetAccess = private)
        dataset % cdm instance
    end

    properties (Dependent = true)
        name            % The string variable name that this object represents
        axes            % the coordinate variables associated with the object
        attributes      % The attributes associated with the object.
    end
    
    properties (SetAccess = private, GetAccess = protected)
       variable 
    end

    
    methods

        %%
        function obj = cdmvariable(cdmInstance, variableName)
            obj.dataset = cdmInstance;
            obj.variable = obj.dataset.dataset.variable(variableName);
        end
        
        function v = get.attributes(obj)
            v = obj.variable.attributes;
        end

        function v = get.axes(obj)
            v = obj.variable.axes;
        end

        function v = get.name(obj)
            v = obj.variable.name;
        end

        function v = size(obj)
            v = obj.variable.size;
        end

        function v = attribute(obj, key)
            v = obj.variable.attribute(key);
        end

        function v = data(obj, varargin)
            [first, last, stride] = obj.toncindex(varargin);
            v = obj.variable.data(first, last, stride);
        end
 
        function v = grid(obj, varargin)
            [first, last, stride] = obj.toncindex(varargin);
            v = obj.variable.grid(first, last, stride);
        end
        
        function v = struct(obj, varargin) 
            [first, last, stride] = obj.toncindex(varargin);
            v = obj.variable.grid(first, last, stride);
            v.(obj.name) = obj.variable.data(first, last, stride);
        end
        
    end
    
    methods (Access = protected)
        function [first, last, stride] = toncindex(obj, a)
            osize = double(size(obj));
            [~, vc] = size(a);
            if vc == 0   % handle .mdata and mdata() cases
                a = {{':'}};
            end
            [first, last, stride] = indexing(a, osize);
        end
        
    end
    
 
    
end
