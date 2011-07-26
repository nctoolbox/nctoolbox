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
% NCTOOLBOX (http://code.google.com/p/nctoolbox)

classdef cfdataset < ncdataset
    
    properties (SetAccess = private, GetAccess = private)
        ncvariables % 2 by n cell array
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
        % function d = data(obj, variable, first, last, stride)
        %    d = squeeze(data@ncdataset(variable, first, last, stride)) 
        % end
        
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
            sz = size(v);
            if (nargin == 2)
              switch length(sz)
                case 1
                  s = v.grid(:);
                  s.(v.name) = v.data(:);
                case 2
                  s = v.grid(:,:);
                  s.(v.name) = v.data(:,:);
                case 3
                  s = v.grid(:,:,:);
                  s.(v.name) = v.data(:,:,:);
                case 4
                  s = v.grid(:,:,:,:);
                  s.(v.name) = v.data(:,:,:,:);
              end
            else
                
                
                % Fill in missing arguments
                % default stride is 1
                if (nargin < 5)
                    stride = ones(1, length(sz));
                end
                
                % Default last is the end
                if (nargin < 4)
                    last = sz;
                end
                
                switch length(sz)
                  case 1
                    s = v.grid(:);
                    s.(v.name) = v.data(first:stride:last);
                  case 2
                    s = v.grid(:,:);
                    s.(v.name) = v.data(first(1):stride(1):last(1),first(2):stride(2):last(2));
                  case 3
                    s = v.grid(:,:,:);
                    s.(v.name) = v.data(first(1):stride(1):last(1),first(2):stride(2):last(2),first(3):stride(3):last(3));
                  case 4
                    s = v.grid(:,:,:,:);
                    s.(v.name) = v.data(first(1):stride(1):last(1),first(2):stride(2):last(2),first(3):stride(3):last(3),...
                      first(4):stride(4):last(4));
                end
                
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
            sz = size(v);
            
            if (nargin == 2)
              switch length(sz)
                case 1
                 s = v.grid(:);
                case 2
                  s = v.grid(:,:);
                case 3
                  s = v.grid(:,:,:);
                case 4
                  s = v.grid(:,:,:,:);
              end
            else
                
                
                % Fill in missing arguments
                % default stride is 1
                if (nargin < 5)
                    stride = ones(1, length(sz));
                end
                
                % Default last is the end
                if (nargin < 4)
                    last = sz;
                end
                
                switch length(sz)
                  case 1
                    s = v.grid(first:stride:last);
                  case 2
                    s = v.grid(first(1):stride(1):last(1),first(2):stride(2):last(2));
                  case 3
                    s = v.grid(first(1):stride(1):last(1),first(2):stride(2):last(2),first(3):stride(3):last(3));
                  case 4
                    s = v.grid(first(1):stride(1):last(1),first(2):stride(2):last(2),first(3):stride(3):last(3),...
                      first(4):stride(4):last(4));
                end
                
            end
        end
        
        function sn = standard_name(obj, standardName)
            % CFDATASET.STANDARD_NAME - Function to get a list of variable names that correspond to
            % to the input standard_name.
            vars = obj.variables;
            for i = 1:length(vars) % Loop through variables to get standard_name attribute if it exists
                tempsn = obj.attribute(vars{i},'standard_name');
                if ~isempty(tempsn) % Make list of standard names, maybe this is not effcient, but it can't
                                               % be too bad?
%                     hash{i, 1} = vars{i};
                    hash{i} = tempsn;
                end
            end
            matches = strcmp(standardName, hash);
            sn = vars(matches);
        end

    end
    
    
    
end

%%
