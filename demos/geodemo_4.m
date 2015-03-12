function geodemo_4
% GEODEMO_4 Compare horizontal slices from two different CF compliant structured
% grid models (CH3D and ROMS) at a particular time step and depth

%% EFDC
url{2}='http://comt.sura.org/thredds/dodsC/data/comt_1_archive/estuarine_hypoxia/VIMS_EFDC/2004_DO3d';
var{2}='salt';
titl{2}='EFDC';

%% ROMS
url{1}='http://comt.sura.org/thredds/dodsC/data/comt_1_archive/estuarine_hypoxia/VIMS_CBOFS/2004-2005';
var{1}='salt';
titl{1}='CBOFS';

dat=[2004 4 10 0 0 0];  % Apr 4, 2004 00:00 UTC

%% Choose a date,depth, and spatial extent
%dat=[2005 1 10 0 0 0];  % Jan 10, 2005 00:00 UTC
depth=-5;  % horizontal slice 5 m from surface
ax=[  -77.4572  -75.4157   36.7224   39.6242]; %lon/lat range
cax=[0 33];  %color range
lat_mid=38; % for scaling plots


%% Perform analysis without using dataset dependant code
% Access datasets, get 3d field at given time data, interpolate data to a constant z, plot results at z depth

figure;
for i=1:length(url);
  nc{i}=ncgeodataset(url{i});
  jd{i}=nj_time(nc{i},var{i});  % using nj_time
  itime=date_index(jd{i},dat);
  disp(['reading data from ' titl{i} '...'])
  
  % using nj_tslice here, which doesn't allow for subsetting or
  % striding:  you get the whole 3D field at a particular time step:
  [s{i},g{i}]=nj_tslice(nc{i},var{i},itime);
  sz{i}=zsliceg(s{i},g{i}.z,depth);
  a{i}=subplot(1,length(url),i);
  pcolorjw(g{i}.lon,g{i}.lat,double(sz{i}));colorbar
  axis(ax);
  caxis(cax);
  title(sprintf('%s, depth=%f: %s',titl{i},depth,datestr(g{i}.time)));
  set (a{i}, 'DataAspectRatio', [1 cos(lat_mid*pi/180) 1000] );

end
