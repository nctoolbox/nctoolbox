function [p1,p2]=iso_plot(data,g,iso_value,vax);
%ISO_PLOT - Plots an isosurface from 3D data and grid structure 
% Usage:[p1,p2] =iso_plot(data,g,iso_value,vax);
% Inputs:  data = 3D data array 
%             g = structure of geo coordinates for data (g.lon,g.lat,g.z) 
%     iso_value = value of isosurface to plot
%     vax = vertical exaggeration (optional, 200 by default)
% Outputs:   p1 = handle of bottom depth surface object
%            p2 = handle of isosurface patch object
% Example: 
% uri='http://geoport.whoi.edu/thredds/dodsC/usgs/vault0/models/examples/bora_feb.nc'
% [t,g]=nj_tslice(uri,'temp',1);% grab 3d field of 'temp' at time step 1
% [p1,p2]=iso_plot(t,g,14,300); % plot 14 degree isosurface with vert exag = 300
%
% NCTOOLBOX (http://code.google.com/p/nctoolbox)
g.z=squeeze(g.z);
data=squeeze(data);
h=squeeze(min(g.z)); % find deepest z values at each grid point
[nz,ny,nx]=size(data);
if (min(size(g.lon))==1),
   [g.lon,g.lat]=meshgrid(g.lon,g.lat);
end
if (min(size(g.z))==1),
    h=h*ones(size(g.lon));
    g.z=repmat(g.z,[1 ny nx]);
end
% plot deepest depth surface 
p1=surf(g.lon,g.lat,h,'FaceColor',[0.7 0.7 0.7],'EdgeColor','none');
if nargin==4,
    fac=111111./vax;    % 111111 meters/degree
else
    fac=111111./200;    % default vertical exxaggeration = 200
end
% plot isosurface
lon3d=permute(repmat(g.lon,[1 1 nz]),[3 1 2]);
lat3d=permute(repmat(g.lat,[1 1 nz]),[3 1 2]);
data(data==0)=nan;
p2 = patch(isosurface(lon3d,lat3d,g.z,data,iso_value));
set(p2, 'FaceColor', [0 1 1], 'EdgeColor', 'none');
daspect([1 cos(mean(g.lat(isfinite(g.lat))*pi/180)) fac]);
% set the view
view(0,45);

%zoom(1.8);
% set lighting
camlight; lighting phong

