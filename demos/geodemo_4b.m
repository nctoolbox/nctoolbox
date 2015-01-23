function geodemo_4b
%% GEODEMO_4b
% Compare horizontal slices from two different CF compliant structured
% grid models (EFDC and CBOFS) at a particular time step and depth, using
% geosubset to subset the data

%% EFDC
url{2}='http://comt.sura.org/thredds/dodsC/data/comt_1_archive/estuarine_hypoxia/VIMS_EFDC/2004_DO3d';
var{2}='salt';
titl{2}='EFDC';

%% ROMS
url{1}='http://comt.sura.org/thredds/dodsC/data/comt_1_archive/estuarine_hypoxia/VIMS_CBOFS/2004-2005';
var{1}='salt';
titl{1}='CBOFS';

%% Create geosubset object
dat=[2004 4 10 6 0 0];  % Apr 4, 2004 06:00 UTC
depth=-5;  % horizontal slice 5 m from surface
ax=[ -76.5220  -75.7105   36.8248   37.7850]; %lon/lat range
cax=[0 33];  %color range
lat_mid=38; % for scaling plots

s.time=datenum(dat);
s.lon=ax(1:2);
s.lat=ax(3:4);

%% Perform analysis without using dataset dependant code
% Access datasets, subset data, interpolate data to a constant z, plot results at z depth

figure;
for i=1:length(url);
  nc{i}=ncgeodataset(url{i});
  % create a salinity geovariable object.  No data read yet.
  svar{i}=geovariable(nc{i},var{i});
  disp(['reading data from ' titl{i} '...'])
  % using geosubset here, which allows for subsetting and striding, reading
  % multiple time steps, only certain zlevels and more.
  sub{i}=svar{i}.geosubset(s);
  sz{i}=zsliceg(squeeze(sub{i}.data),squeeze(sub{i}.grid.z),depth);
  a{i}=subplot(1,length(url),i);
  pcolorjw(sub{i}.grid.lon,sub{i}.grid.lat,double(sz{i}));colorbar
  axis(ax);
  caxis(cax);
  title(sprintf('%s, depth=%f: %s',titl{i},depth,datestr(sub{1}.grid.time)));
  set (a{i}, 'DataAspectRatio', [1 cos(lat_mid*pi/180) 1000] );
end
