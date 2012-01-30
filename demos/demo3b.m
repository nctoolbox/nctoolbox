% DEMO3B - Chaining files and finding attributes using 'value4key'

echo('on')
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
units = value4key(attr, 'units'); % Retrieve the units value
name = value4key(attr, 'long_name'); % Retrieve the long_name value

%% ---- Plot the data
figure;
plot(t, o2);...
datetick('x');...
ylabel([name ' [' units ']']);...
grid;...
shg
echo('off')