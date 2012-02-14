function [index,distance, twoout]=near_lonlat(x,y,x0,y0,dist);
% NEAR_LONLAT  finds the indices of (lon,lat) that are closest to the point (lon0,lat0).
%        [index,distance]=near_lonlat(lon,lat,lon0,lat0) finds the closest point and
%                                           the distance
%        [index,distance]=near_lonlat(lon,lat,lon0,lat0,dist) finds all points closer than
%                                           the value of dist.

% distance=sqrt((x-x0).^2+(y-y0).^2);
distance = zeros(size(x));
for i = 1:length(x(:,1))
    for j = 1:length(x(1,:))
        distance(i,j) = haversine([y(i,j), x(i,j)], [y0, x0]); %this is actually a real distance calculation now...
    end
end
if (nargin > 4),
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