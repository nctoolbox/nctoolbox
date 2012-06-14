function [index, distance, twoout]=near_lonlat2(x,y,x0,y0,dist);
% NEAR_LONLAT finds the indices of (lon,lat) that are closest to the point (lon0,lat0).
%        [index,distance]=near_lonlat(lon,lat,lon0,lat0) finds the closest point and
%                                     the distance(km)
%        [index,distance]=near_lonlat(lon,lat,lon0,lat0,dist) finds all points closer than
%                                     the value of dist(km)
                                       
% Alexander Crosby 2011
% Rich Signell 2012: removed double loop, speeded up 1000x
[nx,ny]=size(x);
x2=x0*ones(size(x));
y2=y0*ones(size(y));
x3=[x(:) x2(:)].';
y3=[y(:) y2(:)].';
distance=sw_dist(y3(:),x3(:),'km');
distance=reshape(distance(1:2:end),nx,ny);

if nargin > 4,
  index=find(distance<=dist);     %finds points closer than dist
  [row col] = find(distance<=dist);
else
  index=find(distance==min(min(distance)));  % finds closest point
  [row col]=find(distance==min(min(distance)));
  index=index(1);
end
distance=distance(index);
twoout = [row col];
end
