%GEODEMO_3  Using ZSLICEG to create a horizontal section from a 3D field
uri='http://geoport.whoi.edu/thredds/dodsC/examples/bora_feb.nc'
[t,g]=nj_tslice(uri,'temp',1);% grab 3d field of 'temp' at time step 1
tz = zsliceg(t,g.z,-5);  % return temperature slice at 5 m depth
figure
pcolorjw(g.lon,g.lat,double(tz));  % plot temp at 5 m depth
title('Horizontal Section of Temperature at 5 m depth');
colorbar