function result=nj_attget(ncRef,var,attrName)
%NJ_ATTGET - Get global and variable attributes
% Usage:   
%   attr= nj_attget(ncRef);                 Returns only 'global' attributes
%   attr= nj_attget(ncRef,'global',attrName)  returns a variable with the value of a global attribute
%   attr= nj_attget(ncRef,var)                Returns only variable attributes
%   attr= nj_attget(ncRef,var,attrName)       returns a variable with the value for a specific variable attribute
% where,
%   ncRef - Reference to netcdf file. It can be either of two
%           a. local file name or a URL  or
%           b. An 'ncgeodataset' matlab object, which is the reference to already
%              open netcdf file.
%              [ncRef=ncgeodataset(uri)]
%   var  = variable name or string 'global'(to retrieve global attributes)   
%   attrName = global or variable attribute name  (case ignored)
% returns,
%   attr - matlab strucure
%
%Sachin Kumar Bhate (skbhate@ngi.msstate.edu) 
%Rich Signell (rsignell@usgs.gov)  (adapted to NCTOOLBOX)

%initialize
result=[];
isNcRef=0;

if nargin < 1, help(mfilename), return, end

%STATIC vars
GLOBAL='global';

try 
    if (isa(ncRef, 'ncgeodataset')) %check for ncgeodataset Object
        nc = ncRef;
        isNcRef=1;
    else
        % open CF-compliant NetCDF File as a Common Data Model (CDM) "Grid Dataset"
        nc = ncgeodataset(ncRef);         
    end       
   
    if (isa(nc, 'ncgeodataset'))
        switch nargin
          case 1      %get all global attributes 
            result=nc.attributes;     
          case 2    %get all variable attributes           
              result = nc.attributes(var);  
           case 3    % get specified global attribute
             if (strcmp(lower(var),'global'))
                 result=nc.attribute(attrName);                
             else % get specified variable attribute    
                 result = nc.attribute(var,attrName);             
             end        
            otherwise, error('MATLAB:nj_attget:Nargin',...
                            'Incorrect number of arguments'); 
        end   

        %cleanup only if we know ncgeodataset object is for single use.
        if (~isNcRef)
            nc.close();
        end
    else
        disp(sprintf('MATLAB:nj_attget:Unable to create "ncgeodataset" object'));
       return;
    end
    
catch 
    %gets the last error generated 
    err = lasterror();    
    disp(err.message); 
end

   
   
end
