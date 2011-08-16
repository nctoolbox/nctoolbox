function cbar = colorbar_log(my_clim)
%COLORBAR_LOG Apply log10 scaling to pseudocolor axis
% and display colorbar COLORBAR_LOG(V), where V is the
% two element vector [cmin cmax], sets manual, logarithmic
% scaling of pseudocolor for the SURFACE and PATCH
% objects. cmin and cmax should be specified on a LINEAR
% scale, and are assigned to the first and last colors in
% the current colormap. A logarithmic scale is computed,
% then applied, and a colorbar is appended to the current
% axis.
%
% Written by Matthew Crema - 7/2007
% Source = Matlab central:
% Subject: colorbar and log scale
% 
% From: Matthew Crema
% 
% Date: 20 Jul, 2007 17:12:42


% Trick MATLAB by first applying pseudocolor axis
% on a linear scale
caxis(my_clim)

% Create a colorbar with log scale
cbar = colorbar('Yscale', 'log');

% Now change the pseudocolor axis to a log scale.
caxis(log10(my_clim));

end
