% DEMO1 - Basic data access
echo('on')
% STARTING DEMO1 ----------------------------------------------------------
% Basic usage of ncdataset

%% ---- Open your data
% Can be NetCDF, NetCDF on webserver, or OpenDAP
ds = ncdataset('http://geoport.whoi.edu/thredds/dodsC/examples/OS_M1_20081008_TS.nc');
%% ---- You can access a list of the variables available to you
ds.variables

%% Lets fetch time in Matlab's native format
t = ds.time('TIME');

%% ---- Now lets get the data
temp = double(ds.data('TEMP')); % Convert to Double!!
depth = ds.data('DEPTH');

%% ---- Plot the data
figure;
surf(t, depth, temp.');...
view(2);shading interp;...
datetick('x', 2);set(gca, 'YDir', 'reverse');...
grid('on');ch = colorbar;...
set(get(ch, 'YLabel'), 'String', '^oC');...
title('Temperature at M1 Mooring in Monterey Bay')
shg

echo('off') % ENDING DEMO1 ------------------------------------------------