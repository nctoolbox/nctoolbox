% DEMO11 - Simple matlab-style use of ncgeodataset
echo on
url = 'http://geoport.whoi.edu/thredds/dodsC/bathy/gom03_v1_0';
nc = cdm(url);
z = nc.data('topo',500:600,400:500);
zg = nc.grid('topo',500:600,400:500);
pcolorjw(zg.lon,zg.lat,z);
% acknowledge source
title({nc.attribute('title'),nc.location,nc.attribute('acknowledgment')},'interpreter','none')
