% DEMO3 - Chaining files and finding attributes

echo('on')
% Starting DEMO3 ----------------------------------------------------------
% An example of chaining data files and finding data attributes (i.e.
% metdata)

o2_ = 'Oxygen';
time_ = 'esecs';

%% ---- Access the OpenDAP datasets
m2006 = ncdataset('http://dods.mbari.org/cgi-bin/nph-nc/data/ssdsdata/deployments/m1/200610/m1_aanderaaoxy_20070105.nc');
m2007 = ncdataset('http://dods.mbari.org/cgi-bin/nph-nc/data/ssdsdata/deployments/m1/200711/m1_aanderaaoxy_20071106.nc');
m2008 = ncdataset('http://dods.mbari.org/cgi-bin/nph-nc/data/ssdsdata/deployments/m1/200810/m1_aanderaaoxy_20081008.nc');
t = [m2006.time(time_); m2007.time(time_); m2008.time(time_)];
o2 = [m2006.data(o2_); m2007.data(o2_); m2008.data(o2_)];

%% ---- Find the units of Oxygen and label the display
attr = m2008.attributes(o2_);
keys = attr(:, 1);    % Cell array of keys
values = attr(:, 2);  % Cell array of values
units = values{ismember(keys, 'units')};  % Retrieve the units value
name = values{ismember(keys, 'long_name')};  % Retrieve the long_name value

%% ---- Plot the data
figure;
plot(t, o2);...
datetick('x');...
ylabel([name ' [' units ']']);...
title('Dissolved Oxygen at M1 Mooring in Monterey Bay');...
grid
shg
echo('off') % Ending DEMO3 ------------------------------------------------