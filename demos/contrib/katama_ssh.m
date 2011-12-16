% KATAMA_SSH
% grab time series of water levels from Chen's GOM3 NECOFS forecast model archive
% for points north and south of Katama Bay, using NCTOOLBOX
%
% Quick NCTOOLBOX Install:
% 1. Grab and unzip: https://github.com/acrosby/nctoolbox_recent/zipball/master
% 2. Run "setup_nctoolbox"

% OPeNDAP URL for Chen's Archive of NECOFS GOM3 Forecast
url='http://www.smast.umassd.edu:8080/thredds/dodsC/fvcom/archives/necofs_gom3';
nc=ncgeodataset(url);
var='zeta';
zvar=nc.geovariable(var);
start=[2011 6 1 0 0 0];
stop=[2011 9 15 0 0 0];

tdat=zvar.timewindowij(start,stop);
% picked these two spots (north, and south of Katama Bay)
loni=[ -70.4893  -70.4873];
lati=[  41.3973   41.3339];
disp(['Reading time series data from ' url '...'])
dat=nj_tseries(nc,'zeta',loni,lati,'method','nearest',...
      'itime',tdat.index,'ele',1);

save katama_ssh.mat dat
% plot whole time series of water levels
figure(1);
set(gcf,'pos',[40 100 900 400]);
plot(dat.time,dat.vals);datetick
legend('north','south');
title('Water levels from NECOFS GOM3')
grid;set(gcf,'color','white');

% zoom into Irene
figure(2)
set(gcf,'pos',[40 150 900 400]);
ii=date_index(dat.time,[2011 8 23 0 0 0],[2011 9 3 0 0 0]);
plot(dat.time(ii),dat.vals(ii,:));datetick
legend('north','south');
title('Water levels from NECOFS GOM3')
grid;set(gcf,'color','white');

%% plot the model bathymetry and locations
figure(3);
lon = nc.data(zvar.getlonname);
lat = nc.data(zvar.getlatname);
depth=nc{'h'}(:);
gvar=zvar.attribute('mesh'); % find mesh variable
grd=nc{gvar}(:); % get the mesh (connectivity array)
[m,n]=size(grd);
if m==3,grd=grd.';elseif n~=3;disp('Error:triangles only');return;end
%zeta=nc{var}(tdat.index(1),:);
%%
trisurf(grd,lon,lat,zeros(size(depth)),-depth);view(2);...
  shading interp;colorbar;dasp(lat(1));
axis([-70.55 -70.40 41.30  41.44]);
caxis([-20 0])
line(loni,lati,'marker','x','linestyle','none','color','black','markersize',14)
title('Katama Bay area bathymetry in NECOFS GOM3 forecast model (m)')
set(gcf,'color','white');
set(gca,'tickdir','out');
