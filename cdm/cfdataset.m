% CFDATASET  Provide access to CF/COARDS convention datasets accessable by the
% NetCDF 4 API. CFDATASET is a subclass of NCDATASET
%
% Use as:
%   ds = cfdataset(dataref)
%
% Arguments:
%   dataref = A reference to a ncdataset that can be accessed by the NetCDF 4
%       API. This includes local netcdf files, netcdf files on web servers
%       and OpenDAP URLs
%
% Return:
%   An instance of a cfdataset class
%   !!! See ncdataset for additional information !!!
%
% Methods:
%   cfdataset.grid - Retrieve coordinate data for the variable.
%   cfdataset.struct - Retrieve data and coordinate data for the variable
%   !!! Refer to ncdataset for additional methods !!!
%
% For more information on the methods use help. For example:
%   >> help cfdataset.struct
%
% Example:
%   ds = cfdataset('http://dods.mbari.org/cgi-bin/nph-nc/data/ssdsdata/deployments/m1/200810/m1_metsys_20081008_original.nc')
%   ga = ds.attributes;       % Global Attributes
%   sv = 'SonicVelocity';     % A variable that we're interested in.
%   d = ds.data(sv);          % Data for the SonicVelocity variable
%   svAx = ds.axes(sv);       % Coordinate Variable names for the SonicVelocity variable
%   svAt = ds.attributes(sv); % Attributes for SonicVelocity
%
% See also NCDATASET

% Brian Schlining
% 2009-10-21
% Alexander Crosby 2010, 2011

classdef cfdataset < ncdataset
    
    properties (SetAccess = private, GetAccess = private)
        ncvariables
    end
    
    methods
        
        %%
        function obj = cfdataset(url)
            % CFDATASET  Constructor. Instantiates a NetcdfDataset pointing to the
            % datasource specified by 'url' and uses that as the underlying
            % dataaccess API. When instantiated, the names of all variables
            % are fetched and stored in the 'variables' property. This can be
            % use to open local files, files stored on an HTTP server and
            % OpenDAP URLs.
            obj = obj@ncdataset(url);
            
            % Java hashTable doens't store matlab objects SO we'll use
            % poor man's hash of a n x 2 cell array
            obj.ncvariables = cell(length(obj.variables), 2);
            for i = 1:length(obj.variables)
                obj.ncvariables{i, 1} = obj.variables{i};
            end
        end
        
        %%
        function v = variable(obj, variableName)
            % CFDATASET.VARIABLE Returns an ncvariable object that provides
            % advanced access to the data contained within that variable.
            %
            % Usage:
            %    v = cfdataset.variable(variableName)
            %
            % Arguments:
            %    variableName = A string name of the variable you want to
            %    retrieve. you can use cfdataset.variables to obtain a list
            %    of all variables available.
            %
            % Returns:
            %
            %    v = an instance of ncvariable
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
                    axesVariableNames = obj.axes(variableName);
                    if ~isempty(axesVariableNames)
                        for i = 1:length(axesVariableNames)
                            if isempty(axesVariableNames{i})
                                axesVariableNames = {};
                                break;
                            end
                        end
                    end
                    
                end
                
                
                
                v = ncvariable(obj, variableName, axesVariableNames);
                if ~isempty(v)
                    for i = 1:length(obj.variables)
                        if strcmp(obj.ncvariables{i, 1}, variableName)
                            obj.ncvariables{i, 2} = v;
                            break;
                        end
                    end
                end
                
            end
            
            
        end
        
        %%
        function s = struct(obj, variableName, first, last, stride)
            % CFDATASET.STRUCT Retrieve all or a subset of the data for the
            % given variable. The data is returned as a structure containing a
            % variable for the data as well as for each dimension of the
            % data.
            %
            % Usage:
            %   d = cfdataset.struct(variableName)
            %   d = cfdataset.struct(variableName, first)
            %   d = cfdataset.struct(variableName, first, last)
            %   d = cfdataset.struct(variableName, first, last, stride)
            %
            %   If no arguments are provided all the data is returned for the
            %   given variable.
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
            %   The data is returned as a structure containing the actual data for the variable
            %   of interest as well as each coordinate variable
            %
            % Example:
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
            
            if (nargin == 2)
                s = v.data;
            else
                sz = size(v);
                
                % Fill in missing arguments
                % default stride is 1
                if (nargin < 5)
                    stride = ones(1, length(sz));
                end
                
                % Default last is the end
                if (nargin < 4)
                    last = sz;
                end
                
                s = v.data(first, last, stride);
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
            %   t = ds.grid('TEMP');
            %
            %
            % Note:
            %   The above example is also equivalent to:
            %   ds = cfdataset('http://dods.mbari.org/cgi-bin/nph-nc/data/ssdsdata/deployments/m1/200810/OS_M1_20081008_TS.nc');
            %   v = ds.variable('TEMP')
            %   t = v.grid
            %   The reason for providing the 'grid' alternative syntax is
            %   that it may be more familiar to some users.
            v = obj.variable(variableName);
            
            if (nargin == 2)
                s = v.grid;
            else
                sz = size(v);
                
                % Fill in missing arguments
                % default stride is 1
                if (nargin < 5)
                    stride = ones(1, length(sz));
                end
                
                % Default last is the end
                if (nargin < 4)
                    last = sz;
                end
                
                s = v.grid(first, last, stride);
            end
        end
        
        %% TODO Commented this out as:
        %  1) gettime makes flawed assumptions on what is a time variable.
        %  2) There may be more than 1 time variable
        %
%         function t = gettime(obj, varargin)
%             t_converted = obj.dataset.time('time'); % this is a bad assumption
%             if nargin > 1
%                 t_index1 = t_converted > varargin{1};
%                 t_index2 = t_converted < varargin{2};
%                 t_index = find(t_index1==t_index2);
%                 t = t_converted(t_index);
%             else
%                 t = t_converted;
%             end
%         end
        
        %%
        % CFDATASET.NUMEL Overridden function required for supporting
        % SUBSREF
        function result = numel(varargin)
            % cfdataset/numel -- Overloaded NUMEL.
            result = 1;
        end
        
        %%
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
                                %                             case '{}' % 5his logic isnt correct
                                %                                 g = substruct('.','variable','()',...
                                %                                   s(2).subs);
                                %                                 d = builtin('subsref',obj,g);
                                %                                 B = builtin('subsref',d,s(3));
                                % %                                 B = d.data(s(3).subs);
                                
                        end
                    else
                        error(['NCTOOLBOX:' mfilename ':subsref'], ...
                            'Unexpected reference notation or method call')
                    end
                    
                 case '()'
                     error('cfdataset:subsref',...
                            'Call with "()" as first type unsupported at this time')
%                     v = obj.variable(char(s(1).subs));
%                     if length(s) == 2
%                         B = v.data(s(2).subs);
%                     elseif length(s) == 3
%                         switch s(3).subs
%                             case 'data'
%                                 B = v.data(s(2).subs);
%                             case 'grid'
%                                 B = v.grid(s(2).subs);
%                         end
%                     else
%                         B = v;
%                     end
                    
                 case '{}'
                    error('cfdataset:subsref',...
                             'Call with "{}" as first type unsupported at this time')       
%                     v = obj.variable(s(1).subs);
%                     if length(s) == 1
%                         B = v;
%                     elseif length(s) == 2
%                         B = v.data(s(2).subs{:});
%                     elseif length(s) == 3
%                         switch s(3).subs
%                             case 'data'
%                                 B = v.data(s(2).subs);
%                             case 'grid'
%                                 B = v.grid(s(2).subs);
%                         end
%                     else
%                         B = v;
%                     end
                    
            end
        end
    end
    
    
    
end

%%
