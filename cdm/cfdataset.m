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
%
%   See ncdataset for additional information.
%
% Methods:
%   Refer to ncdataset for additional methods.
%   cfdataset.variable - Provides advanced access to data for a given variable
%
% For more information on the methods use help. For example:
%   >> help cfdataset.variable
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
                        end
                    end
                end
                
            end
          

        end
        
    end
    
end