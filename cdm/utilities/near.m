function [index,distance]=near(x,x0,n);
% NEAR  finds the indices of x that are closest to the point x0.
%function [index,distance]=near(x,x0,[n]);
%     x is an array, x0 is a point, n is the number of closest points to get
%     (in order of increasing distance).  Distance is the abs(x-x0)
if nargin==2
     n=1;
end
[distance,index]=sort(abs(x-x0));
distance=distance(1:n);
index=index(1:n);

end