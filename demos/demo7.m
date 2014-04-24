% DEMO7 - cfdataset subsetting and axes

% Brian Schlining
% 2009-10-22

echo('on')
% Starting DEMO7 ----------------------------------------------------------
% Demonstration of using cfdataset to subset a variable and it's coordinate axes

%% Open the remote dataset
ds = cfdataset('http://geoport.whoi.edu/thredds/dodsC/examples/OS_M1_20081008_TS.nc');

%% Grab the variable of interest. No data is being read yet.
v = ds.variable('TEMP');

%% Grab a subset of the data. Data is now being pulled across the network
t = v.grid([1 1 1 1], [100 5 1 1], [1 1 1 1]);
t.TEMP = v.data([1 1 1 1], [100 5 1 1], [1 1 1 1]);

%% View the variable name in the returned structure
v.name                  

%% View the coordinate variable names
v.axes        

%% View the data structure.
t                 

%% Make a pretty plot
figure;
plot(ds.time('TIME', t.TIME), t.TEMP) 
cs = cellstr(num2str(t.DEPTH));
for i = 1:length(cs); cs{i} = [cs{i} ' m'];end
legend(cs)
title({ds.attribute('title'),ds.location,['Temperature at ' num2str(t.LATITUDE) 'N, ' num2str(t.LONGITUDE) 'E']},'interpreter','none');
at = ds.attributes('TEMP');
units = value4key(at, 'units');
ylabel(char(units));
xlabel('date');
datetick('x');
grid('on');
shg
echo('off') % Ending DEMO7 ------------------------------------------------
