function [v]=ncload(url,var)
%NCLOAD loads data from all dataset variables (or a single variable) 
% into a structure or directly into the workspace
% Usage: [a]=ncgeodataset(url,[var]);
% Inputs: url = name of dataset (local NetCDF file, Grib, remote OPeNDAP, etc)
%        [var] = optional variable name.  
%                If supplied, loads only that variable
%                If not supplied, loads all variables
% Output: [v] = optional output structure.  
%                If supplied, loads variables into structure
%                If not supplies, loads variables into workspace

% NCTOOLBOX (http://code.google.com/p/nctoolbox)
nc=ncgeodataset(url);
if nargin==1,
  for i=1:length(nc.variables)
    vstr=char(nc.variables(i));
    if nargout==0,
      assignin('base',vstr,nc.data(vstr));
    else
      eval(sprintf('v.%s=nc.data(''%s'');',vstr,vstr));
    end
  end
else
  if nargout==0,
    assignin('base',var,nc.data(var));
  else
    v=nc.data(var);
  end
end
