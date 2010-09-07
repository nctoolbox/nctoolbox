% DEMO1 - Basic data access
echo('on')
% STARTING DEMO1 ----------------------------------------------------------
% Basic usage of ncdataset

%% ---- Open your data
% Can be NetCDF, NetCDF on webserver, or OpenDAP
ds = ncdataset('http://dods.mbari.org/cgi-bin/nph-nc/data/ssdsdata/deployments/m1/200810/OS_M1_20081008_TS.nc');

%% ---- You can access a list of the variables available to you
ds.variables

%% Lets fetch time in Matlab's native format
t = ds.time('TIME');

%% ---- Now lets get the data
temp = double(ds.data('TEMP')); % Convert to Double!!
depth = ds.data('DEPTH');

%% ---- Plot the data
surf(t, depth, temp')
datetick('x', 2);
shading interp
set(gca, 'YDir', 'reverse')
view(2)
grid('on')
ch = colorbar;
set(get(ch, 'YLabel'), 'String', '^oC')
title('Temperature at M1 Mooring in Monterey Bay')

echo('off') % ENDING DEMO2 ------------------------------------------------