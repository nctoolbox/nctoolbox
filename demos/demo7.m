% DEMO7

% Brian Schlining
% 2009-10-22

echo('on')

% Open the remote dataset
ds = cfdataset('http://dods.mbari.org/cgi-bin/nph-nc/data/ssdsdata/deployments/m1/200810/OS_M1_20081008_TS.nc');

% Grab the variable of interest. No data is being read yet.
v = ds.variable('TEMP');

% Grab a subset of the data. Data is now being pulled across the network
t = v.data([1 1 1 1], [100 5 1 1]);

t.metadata.variable    % View the variable name in the returned structure
t.metadata.coordinates  % View the coordinate variable names
t.data                 % View the data variables.

plot(ds.time('TIME', t.data.TIME), t.data.TEMP) 
cs = cellstr(num2str(t.data.DEPTH));
for i = 1:length(cs); cs{i} = [cs{i} ' m'];end
legend(cs)
title(['Temperature at ' num2str(t.data.LATITUDE) 'N, ' num2str(t.data.LONGITUDE) 'E']);

at = ds.attributes('TEMP');
units = value4key(at, 'units');
ylabel(char(units));
xlabel('date');
datetick('x');
grid('on')

echo('off')
