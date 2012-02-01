% GEODEMO_6 
% Demonstrates extracting a layer of velocity from a C-GRID model
% like ROMS
url='http://geoport.whoi.edu/thredds/dodsC/examples/bora_feb.nc';
hname='h'; uname='u'; vname='v'; aname='angle';
nc=ncgeodataset(url);
itime=3; % 3rd time step
klev=-1; % last (top) layer
%Whole domain:
[ U,g ] = cgrid_uv2rho(nc,uname,vname,hname,aname,itime,klev);
figure;
pcolorjw(g.lon,g.lat,abs(U));colorbar;dasp(44);
arrows(g.lon(1:2:end,1:2:end),g.lat(1:2:end,1:2:end),...
  U(1:2:end,1:2:end),0.08,'black');
title(datestr(g.time));dasp(44);
%Subset to jj,ii range:
figure;
[ U,g ] = cgrid_uv2rho(nc,uname,vname,hname,aname,itime,klev,1:58,1:70);
pcolorjw(g.lon,g.lat,abs(U));colorbar;dasp(44);
arrows(g.lon(1:2:end,1:2:end),g.lat(1:2:end,1:2:end),...
  U(1:2:end,1:2:end),0.08,'black');
title(datestr(g.time));dasp(44);
