function zi = griddata_lite(x,y,z,xi,yi,method)
% $Id: griddata_lite.m 358 2008-04-07 14:15:03Z zhang $
% Stripped down version of GRIDDATA that skips a bunch of data checks and 
% saves the delaunay triangulation for use on subsequent calls with the
% exact same set of data coordinates
%
% Invoke griddata_lite the same way as griddata
% e.g. zi = griddata_lite(x,y,z,xi,yi,method)
% 
% Note that only one output is allowed (the interpolated data) and that 
% 'method' must be either 'linear' (the default) or 'nearest'.
%
% Global variable TRI is declared, and if it is empty then the Delaunay
% triangulation is performed as it would be by griddata. The triangulation
% information is saved in global variable TRI (actually a structure) for use 
% on a subsequent call to griddata_lite for input data at the same locations. 
%
% IMPORTANT: If the input data locations alter (e.g. by switching between
% ROMS Arakawa C-grid locations) then the Delaunay triangulation has to be
% redone. This is forced by clearing the global variable TRI and redeclaring 
% it (as empty):
%
%    clear global TRI (be sure to include option global)
%
%
% See help griddata for more info on the method
%
% John Wilkin, October 2002
%              Modified January 2008 to add the 'nearest' option

if nargin < 6
  method = 'linear';
end
switch method
  case { 'linear','nearest'}
  otherwise
    error('only methods LINEAR and NEARSET are supported')
end

global TRI

%%% removed input argument checks %%%

% Need x,y and z to be column vectors
sz = prod(size(x));
x = reshape(x,sz,1);
y = reshape(y,sz,1);
z = reshape(z,sz,1);
%the following added by acrosby to fix for NaN in coordinates
tempindexx = find(isnan(x));
 x(tempindexx)=-[];
 y(tempindexx)=[];
 z(tempindexx)=[];

%%% disabled duplicate data check %%%

% Code lifted from griddata linear and nearest options

%%% function zi = linear(x,y,z,xi,yi)
%%% LINEAR Triangle-based linear interpolation
%   Reference: David F. Watson, "Contouring: A guide
%   to the analysis and display of spacial data", Pergamon, 1994.

siz = size(xi);
xi = xi(:); yi = yi(:); % Treat these as columns
x = x(:); y = y(:); % Treat these as columns

if isempty(TRI)
  disp('      griddata_lite is triangulating data ...')
  % triangularize the data
  tri = delaunayn([x y]);
  if isempty(tri),
    warning('Data cannot be triangulated.');
    zi = repmat(NaN,size(xi));
    return
  end
  switch method
    case 'linear'
      % Find the nearest triangle (t)
      t = tsearch(x,y,tri,xi,yi);
      % Only keep the relevant triangles.
      out = find(isnan(t));
      if ~isempty(out), t(out) = ones(size(out)); end
      tri = tri(t,:);
    case 'nearest'   
      t = [];
      out = [];
  end
  % tri, t, and out are entered into global variable TRI
  TRI.tri = tri;
  TRI.t = t;
  TRI.out = out;
  disp('      ... triangulation complete')
else
  % use the triangulation held in global from last
  % call to griddata_lite
  tri = TRI.tri;
  t = TRI.t;
  out = TRI.out;
end

switch method
  
  case 'linear'
    % Compute Barycentric coordinates (w).  P. 78 in Watson.
    del = (x(tri(:,2))-x(tri(:,1))) .* (y(tri(:,3))-y(tri(:,1))) - ...
      (x(tri(:,3))-x(tri(:,1))) .* (y(tri(:,2))-y(tri(:,1)));
    w(:,3) = ((x(tri(:,1))-xi).*(y(tri(:,2))-yi) - ...
      (x(tri(:,2))-xi).*(y(tri(:,1))-yi)) ./ del;
    w(:,2) = ((x(tri(:,3))-xi).*(y(tri(:,1))-yi) - ...
      (x(tri(:,1))-xi).*(y(tri(:,3))-yi)) ./ del;
    w(:,1) = ((x(tri(:,2))-xi).*(y(tri(:,3))-yi) - ...
      (x(tri(:,3))-xi).*(y(tri(:,2))-yi)) ./ del;
    w(out,:) = zeros(length(out),3);
    
    z = z(:).'; % Treat z as a row so that code below involving
    % z(tri) works even when tri is 1-by-3.
    zi = sum(z(tri) .* w,2);
    zi = reshape(zi,siz);
    if ~isempty(out)
      zi(out) = NaN; 
    end
    
  case 'nearest'
    z = z(:);
    % Find the nearest vertex
    k = dsearch(x,y,tri,xi,yi);
    zi = k;
    d = find(isfinite(k));
    zi(d) = z(k(d));
    zi = reshape(zi,siz);

end
