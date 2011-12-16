function h = pclegend(z,geometry,cfig,label)
% pclegend:  Plot color bar legend anywhere on a plot.
%
% USAGE:
%   h = pclegend( range, geometry, fig_window, label);
%
% PARAMETERS:
%   h:  axis handle of color bar
%
%	PCLEGEND(Z) produces a color bar legend with axis Z.
%	The matrix Z is used to determine the minimum and maximum axis 
%	labels for the color bar. The bar will be horizontal/vertical
%       depending on whether # rows </> # columns of Z.  It is fairly easy
%       to use this if Z = [min_val max_val]'
%
%	PCLEGEND(Z,GEOM) produces the same legend, but at a specific 
%       location on the figure.  The GEOM vector is specified with 
%       [ x y width height ], in normalized coordinates.
%       
%	PCLEGEND(Z,GEOM,FIG,LABEL) does the same with a y axis label and a
%	particular figure window specified.  
 
%       All parameters are optional.  
%
%	Clay M. Thompson  5-28-91
%	Copyright (c) 1991 by the MathWorks, Inc.
%
%       Fudged by me (RP) 9/jan/92 for beta 2, and again on 20/Mar/92 for 
%       beta 3.
%
%       John Evans, 10-26-95
%       Changed to allow for positioning and specified figure windows.

if nargin<1
  z = [0:16]/16; 
end

if nargin < 2
  geometry = [ 0.195 0.15 0.7 0.05 ]; 
end

if nargin < 3
  label = [];
end

if nargin < 4
  cfig = gcf;
end

if nargin>=1
  if isempty(z) 
    z = [0:16]/16; 
  end 
end


[m,n] = size(z);

  zmin = min(z(z~=NaN));
  zmax = max(z(z~=NaN));
  z = zmin + [0:255]*(zmax-zmin)/255;   % Generate 256 equally spaced points

if (n>m),
   clegend_h=axes('position', geometry);
   pcolor(z,[0 10]',[z;z]);
   clegend_h=gca;
   caxis([min(z) max(z)]); 
   shading('flat');
   if nargin>=4, title(label), end
%   if nargin>=4, xlabel(label), end
   set(gca,'tickdir','out','YTick',[]);


else
   clegend_h=axes('position', geometry);
   pcolor([0 10],z,[z;z]');
   clegend_h=gca;
   caxis([min(z) max(z)]); 
   shading('flat');

   if nargin>=4, ylabel(label), end

   set(gca,'tickdir','out','XTick',[]);

%   main_plot_handle = axes('position', [.12 .12 .85-(1-width) .8]);
end;

h = clegend_h;

return;
