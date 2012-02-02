% GEODEMO_6
% Demonstrates extracting a layer of velocity from a C-GRID model
% like ROMS
url = 'http://geoport.whoi.edu/thredds/dodsC/examples/bora_feb.nc';

hname = 'h';
uname = 'u';
vname = 'v';
aname = 'angle';

nc = ncgeodataset(url);

uvar = nc.geovariable(uname);
vvar = nc.geovariable(vname);
hvar = nc.geovariable(hname);
avar = nc.geovariable(aname);

itime = 3; % 3rd time step
klev = -1; % last (top) layer

%Whole domain:
% [ U, g ] = cgrid_uv2rho(nc, uname, vname, hname, aname, itime, klev);
Uobj = hvar.getvectors(uvar, vvar, avar);
g = Uobj.grid(itime, klev, :, :);

% Initialize Plot / Plot Vectors
figure;
pcolorjw(g.lon, g.lat, Uobj.magnitude(itime, klev, :, :));
colorbar;
arrows(g.lon(1:end,1:end), g.lat(1:end,1:end),...
    Uobj.vectors(itime, klev, 1:end,1:end), 0.08, 'black');
title(datestr(g.time));
dasp(44);

%Subset to jj,ii range: New Plot / New Vectors
figure;
% [ U, g ] = cgrid_uv2rho(nc, uname, vname, hname, aname, itime, klev, 1:58, 1:70);
pcolorjw(g.lon(1:58, 1:70), g.lat(1:58, 1:70), ...
    Uobj.magnitude(itime, klev, 1:58, 1:70));
colorbar; 
arrows(g.lon(1:2:58,1:2:70), g.lat(1:2:58,1:2:70),...
    Uobj.vectors(itime, klev, 1:2:58,1:2:70), 0.08, 'black');
title(datestr(g.time));
dasp(44);
