function [h,xx,yy] = pcolorjw(x,y,c)
% $Id: pcolorjw.m 358 2008-04-07 14:15:03Z zhang $
%PCOLORJW
%       PCOLORJW(X,Y,C) is a modified version of PCOLOR that expands the 
%       dimension the inputs X, Y and C so that the checkerboard will show
%       all entries in C.  The default PCOLOR clips the last row and column
%       of C unless the shading('interp') is set, and the color of each 
%       square is determined by the value of C at the squares lower left
%       corner.  
%
%       PCOLORJW recomputes X,Y at the mid points of the input, and pads
%       C so that the effect is to shift each square one half space to the
%       lower left.  The perimeter squares are only half the width/height
%       of all the others as can be seen if one sets shading('faceted')
%
%       John Wilkin 7-May-93
%   
%PCOLOR Pseudocolor (checkerboard) plot.
%       PCOLOR(C) is a pseudocolor plot of matrix C.  The screen is broken
%       into a rectangular "checkerboard" the same size as C.  The value of
%       each element of C specifies the color in its corresponding rectangle.
%       The smallest and largest elements of C are assigned the first and
%       last colors given in the color table; colors for the remainder of the 
%       elements in C are determined by table-lookup within the remainder of 
%       the color table.
%       PCOLOR(X,Y,C), where X and Y are vectors or matrices, makes a
%       pseudocolor plot on the grid defined by X and Y.  X and Y could 
%       define the grid for a "disk", for example.
%       PCOLOR is really a SURF with its view set to directly above.
%       PCOLOR returns a handle to a SURFACE object.
%
%       See also CAXIS, SURF, MESH, IMAGE.

%-------------------------------
%       Additional details:
%
%
%       PCOLOR sets the View property of the SURFACE object to directly 
%       overhead.
%
%       If the NextPlot axis property is REPLACE (HOLD is off), PCOLOR resets 
%       all axis properties, except Position, to their default values
%       and deletes all axis children (line, patch, surf, image, and 
%       text objects).  View is set to [0 90].

%       Copyright (c) 1984-92 by The MathWorks, Inc.

%       J.N. Little 1-5-92

if nargin < 1
        error('Too few input arguments.');
elseif nargin > 4
        error('Too many input arguments.')
end

cax = newplot;
next = lower(get(cax,'NextPlot'));
hold_state = ishold;

if nargin == 2
  error('Must have one or three input arguments.')
end
if nargin == 1
  % error('Use standard PCOLOR')
  c = x;
  [x,y] = meshgrid(1:size(c,2),1:size(c,1));
end
if min([size(x) size(y)]) == 1
  % convert vector inputs to matrices
  [x,y] = meshgrid(x,y);
end
if length(size(c))>2
  % squeeze out the singleton dimension
  c = squeeze(c);
end

% Changes to nc_varget (from snctools) in versions after January 2009 (?) return
% data of the same TYPE as in the actual netcdf file, e.g. FLOAT or DOUBLE or INT.
% Previously data were always converted to DOUBLE. 
% PCOLOR requires DOUBLE type data, so we trap this possibility internaly here
% by simply forcing a type conversion (even if it's not required).
x = double(x);
y = double(y);
c = double(c);

% At this point, x,y and c should all have the same dimension, but
% sometimes c is transposed because of user dyslexia about
% row-column indexing
%
% Trap the case that c needs to be transposed, warn the user, and
% proceed.
if ~all([size(x)-size(c)]==0)
  % size don't match
  if all([size(x)-size(c')]==0)
    % then transposing fixes things
    c = c';
    warning([ 'Input C is being transposed to match X,Y'])
  end
end

[m n] = size(x);
x = [ x(:,1)  0.5*(x(:,1:n-1) + x(:,2:n))  x(:,n)];
y = [ y(:,1)  0.5*(y(:,1:n-1) + y(:,2:n))  y(:,n)];
x = [ x(1,:); 0.5*(x(1:m-1,:) + x(2:m,:)); x(m,:)];
y = [ y(1,:); 0.5*(y(1:m-1,:) + y(2:m,:)); y(m,:)];
c = [ c  NaN*ones(m,1)];
c = [ c; NaN*ones(1,n+1)];
hh = surface(x,y,zeros(size(c)),c);
lims = [min(min(x)) max(max(x)) min(min(y)) max(max(y))];

set(hh,'edgecolor','none','facecolor','flat')

if ~hold_state
  set(cax,'View',[0 90]);
  set(cax,'Box','on');
  axis(lims);
end
if nargout > 0
  h = hh;
  xx = x;
  yy = y;
end
