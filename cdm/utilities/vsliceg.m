function varargout = vsliceg(njData,njGrid,cx,cy)
% VSLICEG:  create a vertical slice along curve
%
% [X,Z,VDATA] = vsliceg(DATA,GRD,CX,CY) creates a vertical slice along
% the coordinates COORDS.  X and Z define the horizontal and vertical
% positions.  It is assumed that DATA and GRD have been returned by
% nj_tslice.
%
% COORDS should be an mx2 array of coordinates contained in GRD.  The first
% column is longitude and the second column is latitude.
%
% VDATA = vsliceg(DATA, GRD, COORDSX, COORDSY) returns just the variable data.
%
% Example:
%    url='http://geoport.whoi.edu/thredds/dodsC/usgs/vault0/models/examples/bora_feb.nc';
%    [data,grd] = nj_tslice(url,'temp',1);
%    pcolorjw(grd.lon,grd.lat,double(data(end,:,:)));
%    [cx,cy] = ginput(10); %click on 10 points
%    [x,y,vdata] = vsliceg(data,grd,cx,cy);
%    pcolorjw(x,y,vdata);
%    shading flat
%
%

error(nargchk(4,4,nargin,'struct'));
coords=[cx(:) cy(:)];

num_levels = size(njGrid.z,1);

% Force the number of slice points to be at least 100.
numpts = size(coords,1);
if numpts < 100
    t = linspace(0,1,numpts);
    ft = linspace(0,1,100);
    slice_x = spline(t,coords(:,1),ft);
    slice_y = spline(t,coords(:,2),ft);
else
    slice_x = coords(:,1);
    slice_y = coords(:,2);
end

% need to expand 1D lon,lat to 2D lon,lat
% need to expand 1D z to 3D z
[nz,ny,nx]=size(njData);
if (min(size(njGrid.lon))==1),
    [njGrid.lon,njGrid.lat]=meshgrid(njGrid.lon,njGrid.lat);
end
if (min(size(njGrid.z))==1),
    njGrid.z=repmat(njGrid.z,[1 ny nx]);
end


% We need to find out if any of the slice points are outside of
% the grid.  Data for these points shall be NaN.
njGrid.lon=double(njGrid.lon);
njGrid.lat=double(njGrid.lat);

[nr,nc] = size(njGrid.lon);
boundary_lon = [njGrid.lon(1:nr-1,1); njGrid.lon(nr,1:nc-1)'; njGrid.lon(nr:-1:2,nc); njGrid.lon(1,nc:-1:2)'];
boundary_lat = [njGrid.lat(1:nr-1,1); njGrid.lat(nr,1:nc-1)'; njGrid.lat(nr:-1:2,nc); njGrid.lat(1,nc:-1:2)'];
vslice_pt_inside_njGrid = inpolygon ( slice_x, slice_y, boundary_lon, boundary_lat );

% Form the horizontal distance along the slice.
r = geodist_wr ( slice_x, slice_y );
x = [0; cumsum(r)];

% Convert to kilometers
x = x/1000;
x = repmat ( x', num_levels, 1 );
transect_dist = x(1,:);


% Interpolate the data using GRIDDATA_LITE

% This method relies upon a global variable called TRI.  We need
% to clear it each time this is called.
clear global TRI;
global TRI; %#ok<NUSED>


% Compute the triangulation exactly once, the first time thru.
%
% 3D field.  Have to interpolate depths as well.
clear global TRI;
global TRI; %#ok<REDEF>
vdata = zeros(num_levels,numel(transect_dist));
vdata(1,:) = griddata_lite ( njGrid.lon, njGrid.lat, squeeze(njData(1,:,:)), slice_x, slice_y );
for zi = 2:num_levels
    vdata(zi,:) = griddata_lite ( njGrid.lon, njGrid.lat, squeeze(njData(zi,:,:)), slice_x, slice_y );
end



% Grid the proper depths.
% This seems to assume that grid points don't vary with depth.
clear global TRI;
global TRI;
vz(1,:) = griddata_lite ( njGrid.lon, njGrid.lat, squeeze(njGrid.z(1,:,:)), slice_x, slice_y );
for zi = 2:num_levels
    vz(zi,:) = griddata_lite ( njGrid.lon, njGrid.lat, squeeze(njGrid.z(zi,:,:)), slice_x, slice_y );
end

vz = flipud(vz);

% i dont know how important this is, but commenting it out fixes a
% problem...
% The last thing to do is to NaN out any of the data points that fell outside the
% grid.
% n = length(transect_dist);
% for j = 1:n
%     if ~vslice_pt_inside_njGrid(j)
%         vdata(:,j) = NaN;
%     end
% end
vdata = flipud ( vdata );




switch(nargout)
    case 1
        varargout{1} = vdata;
    case 3
        varargout{1} = x;
        varargout{2} = vz;
        varargout{3} = vdata;
end



%-------------------------------------------------------------------------------
function r = geodist_wr ( lon, lat )
% GEODIST_WR:  distance between lat/lon coords in meters..
%
% uses the geodetic inverse formula (geodesy by %  A.R. Clarke)
% Actually, this is just a wrapper around geodist.  It agrees
% with what the "distance" routine in the matlab mapping toolbox
% thinks to within maybe 1 percent or so.  Good enough for
% government work.
%
% USAGE:  r = geodist_wr ( lon, lat );
%
% PARAMETERS:
% Input:
%     lon, lat:
%         in degrees
% Output
%     r:
%         distance (range) between points in meters


n = length(lon);
r = zeros(n-1,1);
for j = 1:n-1
    r(j,1) = geodist ( [lat(j) lon(j) lat(j+1) lon(j+1)] );
end

return


function [d,az1,az2]=geodist(input)
%  GEODIST uses the geodetic inverse formula (geodesy by
%  A.R. Clarke)
%  input is a vector in the form (slat slon elat elon) such that
%            slat, slon - station latitude and longitude (degrees)
%            elat, elon - shot latitude and longitude (degrees)
%  output is:  d - distance in meters
%              az1 - azimuth station from shot (degrees)
%              az2 - azimuth shot from station (degrees)
%  adapted from a fortran subroutine used in the GLIMPCE expt.
%
%  D. Hutchinson  1 September 1994






able = 6378.2064;
bake = 6356.5838;
easy = .006768658;

slat=input(1)*pi/180  ;                          %convert to radians
slon=input(2)*pi/180;
elat=input(3)*pi/180;
elon=input(4)*pi/180;

rat=(bake*bake)/(able*able);
cslat=cos(slat);
sslat=sin(slat);
celat=cos(elat);
selat=sin(elat);
cd=cos(slon-elon);
sd=sin(slon-elon);
ens=able/sqrt(1.-easy*sslat*sslat);
ene=able/sqrt(1.-easy*selat*selat);
ene=ene/ens;
a=-cslat*sd;
b= rat*celat*sslat -cslat*selat*cd +easy*celat*selat*ene;
az1= atan2(a,b);
a = celat*sd;
b = rat*cslat*selat -celat*sslat*cd +easy*cslat*sslat/ene;
az2 = atan2(a,b);
fkon = celat*cd*ene - cslat;
a = celat*sd*ene;
b = selat*ene - sslat;
b = b*rat;
fkon = fkon*fkon + a*a +b*b;
fkon = sqrt(fkon);
b = sqrt(easy/(1-easy));
a = b*sslat;
b = b*cslat*cos(az2);
beff = 1 + b*b;
fkor = fkon*beff;
beff = 1/beff;
bach = (a*a - b*b)*beff;
beff = a*b*beff;
c = zeros(1,4);
c(1) = -.1875*beff*bach - (1/3)*beff*beff*beff;
c(2) = .0046875 + .0375*bach + .25*beff*beff;
c(3) = -.125*beff;
c(4) = .04166666667;
d = c(1);
for k = 2,4;
    d = d*fkor + c(k);
end

d = d*fkor*fkor + 1;
d = (fkon*ens*d)*1000;
az1 = az1*180/pi;
az2 = az2*180/pi;





