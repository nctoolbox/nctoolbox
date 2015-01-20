function a = getaxes(netcdf, variableName, useSimple)



end

%%
function a = simpleAxes(netcdf, variableName)
% NCDATASET.AXES  Returns a cell array containing the variable names of
% coordinate axes for the given variable
%
% Use as:
%   ds = ncdataset('http://dods.mbari.org/cgi-bin/nph-nc/data/ssdsdata/deployments/m1/200810/m1_metsys_20081008_original.nc')
%   ax = ds.axes(variableName);
%
% Inputs:
%   variableName = the name of the variable whos axes you want to retrieve
%
% Return:
%   An (n, 1) cell array containing the names (in order) of the variable
%   names representing the coordinate axes for the specified variableName.
v = obj.findvariable(variable);
dims = v.getDimensions();
cv = cell(dims.size(), 1);
for i = 1:dims.size();
    ca = obj.netcdf.findCoordinateAxis(dims.get(i - 1).getName());
    if ~isempty(ca)
        cv{i} = char(ca.getName());
    end
end
end