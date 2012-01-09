%coawst_isosurface_now
% isosurface movie from COAWST (ROMS) forecast in Gulf of Maine
% plot temperature isosurface every 25 hours over the last 90 days

% COAWST (forecast archive + current forecast)
url='http://geoport.whoi.edu/thredds/dodsC/coawst_2_2/fmrc/coawst_2_2_best.ncd';
var='temp';
movie_name='coawst_iso'
isoval=10;  % isosurface_value
nhours=90*24;  % how many hours back from current time
nsub=25;  % subsample in hours
titl='COAWST 10 degree Temperature Isosurface';
vert_exag=200; % vertical exaggertation for plot
lon_range=[-71.4 -63];lat_range=[41 46];
hstride=2; % stride in lon/lat space (e.g. hstride=2 => every other point)
jd=now_utc; % gets the current time in UTC

disp(['Opening dataset ' url '...'])
nc=ncgeodataset(url);
jd=nj_time(nc,var);   % hourly data
ii=[(length(jd)-nhours):nsub:length(jd)];
uvar=nc.geovariable(var);
vidobj=VideoWriter(movie_name); %creates AVI file, test.avi 
open(vidobj);
fig=figure(1);
urlc='http://geoport.whoi.edu/thredds/dodsC/usgs/vault0/data/coast/gom_coast.nc';
ncc=ncgeodataset(urlc);lonc=ncc{'lon'}(:);latc=ncc{'lat'}(:);

for i=1:length(ii)
  disp(['Extracting subsetted data...'])
  iii=ii(i);
  s.t_index=iii;
  s.lon=lon_range;
  s.lat=lat_range;
  s.h_stride=[hstride hstride]
  tic;
  if(i==1)
    us=uvar.geosubset(s);
    data=us.data;
    grid=us.grid;
    ind=us.indices;
    tocs(i)=toc;
  else
    data=nc{var}(iii,ind{2},ind{3},ind{4});
    tocs(i)=toc;
  end
  % call iso_plot (NCTOOLBOX routine)
  iso_plot(data,grid,isoval,vert_exag)
  % coastline
  line(lonc,latc,zeros(size(lonc)),'color','black');
  % set camera just the way we want it (I got this by zooming in manually
  % and then clicking "File=>Generate Code" in the figure window)
  set(gca,...
  'PlotBoxAspectRatio',[1.16666666666667 1.35866104619474 0.675000675000675],...
  'DataAspectRatio',[1 0.736018746397965 555.555],...
  'CameraViewAngle',3.3377073141751,...
  'CameraUpVector',[0 0 1],...
  'CameraTarget',[-65.7780743350319 41.0728486210621 -1230.73456512896],...
  'CameraPosition',[-136.210169005214 74.7377109586742 42276.5651353494]);
  title(sprintf('%s:%s',titl,datestr(jd(iii))),'position',[-72.3 43.9 0]);
  F=getframe(fig);
  %anim_frame(movie_name,jj)
  writeVideo(vidobj,F);
end
close(vidobj);