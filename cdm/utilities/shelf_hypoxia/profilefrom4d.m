function [temp lat lon] = profilefrom4d(temp, out)
%
%
%
%
%
if length(temp.var_struct.grid.lon) ~= length(temp.var_struct.grid.lat);
    [temp.var_struct.grid.lon temp.var_struct.grid.lat] = meshgrid(temp.var_struct.grid.lon, temp.var_struct.grid.lat);
else
end
if length(size(temp.var_struct.grid.z)) > 2
    
    temp.depths = double(squeeze(temp.var_struct.grid.z(:, nearxy(temp.var_struct.grid.lon, temp.var_struct.grid.lat, out(1), out(2)))));
else
    
    temp.depths = temp.var_struct.grid.z;
end

temp.profile = double(squeeze(temp.var_struct.data(1, :, nearxy(temp.var_struct.grid.lon, temp.var_struct.grid.lat, out(1), out(2)))));

lat = temp.var_struct.grid.lat(nearxy(temp.var_struct.grid.lon, temp.var_struct.grid.lat, out(1), out(2)));
lon = temp.var_struct.grid.lon(nearxy(temp.var_struct.grid.lon, temp.var_struct.grid.lat, out(1), out(2)));
end