% NCUVARIABLE Extention of variable object class for unstructured grid/mesh variables.
% Trying out sort of an interface approach, I think I like it - Acrosby
classdef ncuvariable < handle
    
    properties (SetAccess = private)
        dataset % matlab ncugrid class
        variable %netcdf meshvariable
        name
        axes
        size
    end
    
    methods
        
        function obj = ncuvariable(src, variableName)
            % Constructor
            if ischar(src)
                obj.dataset = ncugrid(src);  % src is a string URL/File
            elseif isa(src, 'ncugrid')
                obj.dataset = src;             % src is a ncdataset
            else
                ex = MException('NCTOOLBOX:ncuvariable', 'Invalid dataset was specified');
                ex.throw;
            end
            
            obj.variable = src.netcdfugrid.getMeshVariableByName(variableName);
            
            obj.name = variableName;
            
            coordinates = src.attribute('coordinates', variableName);
            coordinate_vars = regexp(coordinates, ' ', 'split')';
            obj.axes = coordinate_vars;
            
            obj.size = src.size(variableName);
            
        end % constructor end
        
        function s = subset(obj, varName, struct);
            s = obj.dataset.unstructuredLatLonSubset(varName, struct);
        end % subset end
        
        function atts = attributes(obj)
            atts = obj.dataset.attributes(obj.name);
        end % attributes end
        
        function att = attribute(obj, key)
            att = obj.dataset.attribute(key, obj.name);
        end % attribute end
        
        function d = data(obj, first, last, stride)
            d = obj.dataset.data(obj.name, first, last, stride);
        end % data end
        
        function g = grid(obj, first, last, stride)
            g = obj.dataset.grid(obj.name, first, last, stride);
        end %grid end
        
        function gi = grid_interop(obj, first, last, stride)
            gi = obj.dataset.grid_interop(obj.name, first, last, stride);
        end %grid_interop end
        
        %%
        function e = end(obj, k, n)
            n = obj.size;
            e = n(k);
        end % Added to deal with end indexing functionality,
        % otherwise the indexing arugment is ignored.
        
        function sref = subsref(obj,s)
            switch s(1).type
                % Use the built-in subsref for dot notation
                case '.'
                    switch s(1).subs
                        case 'grid_interop'
                            switch length(s)
                                case 1
                                    sref = obj;
                                    %                           obj.temp = size(s(2).subs{1});
                                case 2
                                    nums = obj.size;
                                    [first last stride] = indexing(s(2).subs, double(nums));
                                    sref = obj.grid_interop(first, last, stride);
                            end
                        case 'data'
                            
                            %                             if ~isempty(nums)
                            switch length(s)
                                case 1
                                    sref = obj;
                                    %                                         obj.temp = size(s(2).subs{1});
                                case 2
                                    nums = obj.size;
                                    [first last stride] = indexing(s(2).subs, double(nums));
                                    sref = obj.data(first, last, stride);
                            end
                            
                            %                             else
                            %                                 sref = obj.data;
                            %                                 warning(['NCTOOLBOX:ncgeovariable:subsref'], ...
                            %                                     ['Variable "' name '" has no netcdf dimension associated with it. Errors may result from non CF compliant files.'])
                            %                             end
                        case 'grid'
                            
                            %                             if ~isempty(nums)
                            switch length(s)
                                case 1
                                    sref = obj;
                                    %                                         obj.temp = size(s(2).subs{1});
                                case 2
                                    nums = obj.size;
                                    [first last stride] = indexing(s(2).subs, double(nums));
                                    sref = obj.grid(first, last, stride);
                            end
                            
                            %                             else
                            %                                 warning(['NCTOOLBOX:ncgeovariable:subsref'], ...
                            %                                     ['Variable "' name '" has no netcdf dimension associated with it. Errors may result from non CF compliant files.'])
                            %                             end
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
        end % subsref end
        
    end % methods end
    
end % class end