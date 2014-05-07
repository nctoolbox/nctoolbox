% DEMO12 - Simple matlab-style use of ncgeodataset{varname}() syntax
echo on
url='http://tds.marine.rutgers.edu:8080/thredds/dodsC/projects/wilkin/southchinasea/his/luzonrun003/his_scs2km_003_0027.nc'

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

