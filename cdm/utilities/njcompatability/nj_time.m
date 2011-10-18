function result=nj_time(ncRef,var,timeIndex)
% NJ_TIME - Get DATENUM times for a variable with CF-compliant time
% Usage:
%   result=nj_time(ncRef,var);   - all times
%   result=nj_time(ncRef,var,timeIndex); - time for particular timestep
% where,
%   ncRef - Reference to netcdf file. It can be either of two
%           a. local file name or a URL  or
%           b. An 'ncgeodataset' matlab object, which is the reference to already
%              open netcdf file.
%              [ncRef=ncgeodataset(uri)]
%   var - variable name
%   timeindex - time step
% Returns,
%   result - vector of times in DATENUM format
%
% e.g,
%  Gridded dataset:
%       uri='http://coast-enviro.er.usgs.gov/cgi-bin/nph-dods/models/test/bora_feb.nc';
%       var='temp';
%  Non-Gridded dataset: (But has CF-compliant time coordinate variable)
%       uri='http://www.gri.msstate.edu/rsearch_data/nopp/non-gridded/adcirc/adc.050913.fort.64.nc';
%       var='ubar';
%  
%
% Richard Signell (rsignell@usgs.gov)
% Sachin Kumar Bhate (skbhate@ngi.msstate.edu) 

%initialize
result=[];
isNcRef=0;

if nargin < 2, help(mfilename), return, end

try 
    if (isa(ncRef, 'ncgeodataset')) %check for ncgeodataset Object
        nc = ncRef;
        isNcRef=1;
    else
        % open CF-compliant NetCDF File as a Common Data Model (CDM) "Grid Dataset"
        nc = ncgeodataset(ncRef);         
    end
    
     if (isa(nc, 'ncgeodataset'))
        ncVar = nc.geovariable(var);
        time_name=ncVar.gettimename();
        switch nargin              
          case 2    %variable attributes 
              result = nc.time(time_name);  %all times                
          case 3    % get specific time step
              var_grid = ncVar.grid_interop(timeIndex);
              result = var_grid.time;      
           otherwise, error('MATLAB:nj_time:Nargin',...
                            'Incorrect number of arguments'); 
        end 
    else
        disp(sprintf('MATLAB:nj_time:Unable to create "ncgeodataset" object'));
       return;
    end 
catch 
    %gets the last error generated 
    err = lasterror();    
    disp(err.message); 
end
end



