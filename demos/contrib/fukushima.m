% FUKUSHIIMA get surface currents from Fukushima ocean forecast via OPeNDAP

% NCTOOLBOX:  http://code.google.com/p/nctoolbox/
url='http://edac-dap3.northerngulfinstitute.org/thredds/dodsC/ncom_fukushima_agg/Fukushima_best.ncd';
uvar='water_u';
vvar='water_v';
ksurface=1;
start=[2011 3 20 0 0 0];  % want data closest to 2010-Mar-20 19:45:00Z
stop=[2011 3 20 12 0 0];   %subsample lon,lat?
%isub=1; %no subsampling in lon,lat
isub=2;  %subsample to speed up test
ivec=4; % additional subsampling of vectors
disp(['Opening ' url '...'])
nc=ncgeodataset(url);
disp(['Reading data...'])
[lon,lat]=nj_lonlat(nc,uvar);
lon=lon(1:isub:end);
lat=lat(1:isub:end);
[lon2d,lat2d]=meshgrid(lon(1:ivec:end),lat(1:ivec:end));
jd=nj_time(nc,uvar);
ind=date_index(jd,start,stop);
for i=1:length(ind)
  udata=nc{'water_u'}(ind(i),ksurface,1:isub:end,1:isub:end);
  vdata=nc{'water_v'}(ind(i),ksurface,1:isub:end,1:isub:end);
  w=double(squeeze(complex(udata,vdata)));
  pcolor(lon,lat,abs(w));shading flat;colorbar;dasp(35);
  title(sprintf('NCOM Fukushima forecast surface currents (m/s): %s',...
    datestr(jd(ind(i)))));
  caxis([0 2]);set(gca,'tickdir','out');grid;
  arrows(lon2d,lat2d,w(1:ivec:end,1:ivec:end),0.08,'black');
  shg
end
