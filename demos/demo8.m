% DEMO8

echo('on')
% Starting DEMO8 ----------------------------------------------------------
% Demonstration of subsetting a CF convention dataset

url='http://geoport.whoi.edu/thredds/dodsC/coawst_2_2/fmrc/coawst_2_2_best.ncd'; 
ds = cfdataset(url);

% Grab the variable of interest. No data is being read yet.
% NOTE: You should use ds.struct now instead of ds.variable. See demo9.m
v = ds.variable('temp');
% sz= size(v);

% Grab a subset of the data. Data is now being pulled across the network
t = v.data(end, end, 1:end, 1:end);
g = v.grid(end, end, 1:end, 1:end);

% Make a pretty plot. Note the call to 'squeeze'. This removes
% singleton dimensions.
figure;
surf(g.lon_rho, g.lat_rho, double(squeeze(t)))
shading('interp');
view(2)
axis('equal')

xatt = ds.attributes('lon_rho');
xname = value4key(xatt, 'long_name');
xunits = value4key(xatt, 'units');
xlabel([xname ' [' xunits ']']);

yatt = ds.attributes('lat_rho');
yname = value4key(yatt, 'long_name');
yunits = value4key(yatt, 'units'); 
ylabel([yname ' [' yunits ']']);

zatt = ds.attributes('temp');
zname = value4key(zatt, 'long_name');
zunits = value4key(zatt, 'units'); 
ztime = ds.time('time1', g.time1);
title([zname ' [' zunits '] on ' datestr(ztime)]);

colorbar
shg

echo('off') % Ending DEMO8 ------------------------------------------------