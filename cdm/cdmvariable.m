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
classdef cdmvariable < cfvariable
    
    
    methods
        
        %%
        function obj = cdmvariable(src, variableName, axesVariableNames)
            % CFVARIABLE.CFVARIABLE  Constructor.
            %
            % Use as:
            %    v = cfvariable(src, variableName)
            %    v = cfvariable(src, variableName, axesVariableNames)
            %
            obj = obj@cfvariable(src, variableName, axesVariableNames)
            
            
        end
        
 
 
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
                                        sref = obj.data;
                                    case 2
                                        idx = s(2).subs;
                                        [first, last, stride] = indexing(idx, obj.size);
                                        sref = obj.somedata(1, first, last, stride);
                                end
                                
                            else
                                sref = obj.data;
                                warning('NCTOOLBOX:cfvariable:subsref', ...
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
                                        [first, last, stride] = indexing(idx, obj.size);
                                        sref = obj.grid(first, last, stride);
                                end
                                
                            else
                                warning('NCTOOLBOX:cfvariable:subsref', ...
                                    ['Variable "' name '" has no netcdf dimension associated with it. Errors may result from non CF compliant files.'])
                            end
                        otherwise
                            sref = builtin('subsref',obj,s);
                    end
                case '()'
                    if length(s) < 2
                        % Note that obj.Data is passed to subsref
                        sref = builtin('subsref', obj.data, s);
                    else
                        sref = builtin('subsref', obj, s);
                    end
                    % No support for indexing using '{}'
                case '{}'
                    error('NCTOOLBOX:cfvariable:subsref', ...
                        'Not a supported subscripted reference, "{}" are not permitted to call variable object methods');
            end
        end
        
        
    end
    
 
    
end
