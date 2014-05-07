% DEMO12 - Simple matlab-style use of ncgeodataset{varname}() syntax
echo on
url='http://geoport.whoi.edu/thredds/dodsC/usgs/vault0/models/examples/bora_feb.nc'
nc=ncgeodataset(url)
varname = 'temp';
V=nc.geovariable(varname);
lon=nc{V.getlonname}(:,:); 
lat=nc{V.getlatname}(:,:); 
temp=nc{varname}(1,1,:,:);
pcolorjw(lon,lat,temp);
colorbar
% acknowledge source
title({nc.attribute('title'),nc.location,V.attribute('long_name')},'interpreter','none')

