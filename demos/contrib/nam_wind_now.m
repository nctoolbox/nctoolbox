% NAM_WIND_NOW
% Extract and plot data from NAM 4km/12km/33km data from NCEP NCO NOMADS.
% This plots current wind in the Gulf of Maine, but could easily be 
% modified for different variables, regions, times...
% Info: NOMADS page: http://nomads.ncep.noaa.gov/index.shtml
% NAM OPeNDAP: http://nomads.ncep.noaa.gov:9090/dods/nam
%   then select directory by date "e.g. nam20111120" (after 20111114)
%    - select "nam_conusnest_00z" or "nam_conusnest_06z", for 4km
%    - select "nam_00z" or "nam_06z", for 12 km
%    - select "nam_na_00z" or "nam_na_06z", for 33 km

% Rich Signell (rsignell@usgs.gov)
% NCTOOLBOX (http://code.google.com/p/nctoolbox/)

jd=now_utc; % gets the current time in UTC
g=datevec(jd-6/24); % go for forecast 6 hours before
d=datestr(g,'yyyymmdd');
h=floor(g(4)/6)*6;

% NAM 4km
url=sprintf('http://nomads.ncep.noaa.gov:9090/dods/nam/nam%s/nam_conusnest_%2.2dz',d,h);
lon_range=[-71.5 -63];lat_range=[41 46];  % 4km lon is [-180 180]

% NAM 12km
%url=sprintf('http://nomads.ncep.noaa.gov:9090/dods/nam/nam%s/nam_%2.2dz',str,d,h);
%lon_range=[-71.5 -63];lat_range=[41 46];  % 12km lon is [-180 180]

% NAM 32km
%url=sprintf('http://nomads.ncep.noaa.gov:9090/dods/nam/nam%s/nam_na_%2.2dz',str,d,h);
%lon_range=[-71.5 -63]+360;lat_range=[41 46]; % 32km lon is [0 360]

nc=ncgeodataset(url);
disp(['Opening dataset ' url '...'])

uvar=nc.geovariable('ugrd10m');
vvar=nc.geovariable('vgrd10m');

%s.time=jd; % datenum is simpler, but datevec is more readable
%s.time=[{datevec(jd)}];
s.time=datestr(jd);
s.lon=lon_range;
s.lat=lat_range
disp(['Reading data...'])
us=uvar.geosubset(s);
vs=vvar.geosubset(s);
lon=us.grid.lon;
lat=us.grid.lat;
w=squeeze(complex(us.data,vs.data)); % vector wind as complex
pcolorjw(lon,lat,abs(w));
dasp(44);
colorbar
title(['NCEP NAM: 10m Wind Speed (m/s): '  datestr(us.grid.time) ' UTC'])
set(gca,'tickdir','out');set(gcf,'color','white');
% subsample arrows by isub
[xx,yy]=meshgrid(lon,lat);
isub=2;
g=arrows(xx(1:isub:end,1:isub:end),yy(1:isub:end,1:isub:end),...
  w(1:isub:end,1:isub:end),0.01,'black');
shg

% load & plot a custom Gulf of Maine coastline: 
% this is non-essential, and can be commented out
urlc='http://geoport.whoi.edu/thredds/dodsC/usgs/vault0/data/coast/gom_coast.nc';
ncc=ncgeodataset(urlc);lonc=ncc{'lon'}(:);latc=ncc{'lat'}(:);
line(lonc,latc,'color','black');


