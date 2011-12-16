% Script to access surface currents from Fukushima ocean forecast
% via OPeNDAP using NCTOOLBOX 
% http://code.google.com/p/nctoolbox/

%  This script uses latest version from Mecurial, which contain
%  a new matlab-based indexing not yet in the public release:
%  hg clone https://nctoolbox.googlecode.com/hg/ nctoolbox 
%  Rich Signell (rsignell@usgs.gov)


%WARNING TO MATT: Some datasets may not be in the correct format, if this is the case some errors may occur
% If ncgeodataset/ geovariable are giving you problems try replacing ncgeodataset with 'cfdataset' and nc.geovariable with 'nc.variable'.
% If you do replace the classes you are using some functions like grid_interop won't work, but they have slightly less sophisticated versions like
% nc.grid.


url='http://edac-dap3.northerngulfinstitute.org/thredds/dodsC/ncom_fukushima_agg/Fukushima_best.ncd';
nc=ncgeodataset(url); % Get dataset object

vars = nc.variables % List variable names. nc.attributes will list the dataset attributes and
                             % nc.attributes('variablename') will list variable attributes  
                             
nc.location % Is where the path or link is pointing to incase you forget what dataset you were working with

nc.size('water_u') % Gets the size/dimensions of a variable if you want to look at it (before you get the values themselves, otherwise you can just do the standard size(umatrix) thing)

u=nc.geovariable('water_u'); % Get objects for variables of interest
                                           % u.size will also get the size of the variable
v=nc.geovariable('water_v');

ksurface=1; % Just the vertical level of interest

isub=4;  %optionally subsample lon/lat array for speed/testing

grd=u.grid_interop(:, 1, 1:isub:end, 1:isub:end);  % Returns a matlab structure with values corresponding to coordinates of a variable
dn = grd.time;  % toolbox already converts whatever format time is in to matlab datenum format

%for itime=1:length(dn);
for itime=1:10;
    % This calls a slice of data for every loop. Sometimes you want the entire variable as one chunk, but sometimes it makes sense to do it bit by
    % bit, especially if you hit server limits. (403 Type errors)
    udata=squeeze(u.data(itime, ksurface, 1:isub:end, 1:isub:end)); % Data is 4d but should be 2d for plotting since 2 of the dims are length of 1
    vdata=squeeze(v.data(itime, ksurface, 1:isub:end, 1:isub:end)); % Data is 4d but should be 2d for plotting since 2 of the dims are length of 1
    
    
    disp(sprintf('Step %d/%d',itime,length(dn)))
    
    % This makes data set complex/vector for finding the magnitude with abs() in the next step
    w=double(complex(udata,vdata)); % Sometime model data is a 'single' type instead of 'double' so we have to convert
    
    
    pcolor(grd.lon,grd.lat,abs(w)); % x and y values come from the grid_interop command above
                                                  % We know that u and v are on the same grid so we only get the coordinates once, but that may not
                                                  % always be the case
    
    shading flat;
    colorbar
    
    title(sprintf('NCOM Fukushima forecast surface currents (m/s): %s',...
        datestr(dn(itime))));
    caxis([0 2]);set(gca,'tickdir','out');
    drawnow
end