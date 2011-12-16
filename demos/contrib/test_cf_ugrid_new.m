% Test script to read CF-UGRID convention data from
% ADCIRC, FVCOM and ELCIRC, and SELFE
%
uris={'dods://testbedapps.sura.org/thredds/dodsC/ugrid/TestCases/ADCIRC/adcirc_delt.ncml',...
    'dods://testbedapps.sura.org/thredds/dodsC/ugrid/TestCases/FVCOM/fvcom_delt.ncml',...
    'dods://testbedapps.sura.org/thredds/dodsC/ugrid/TestCases/ELCIRC/elcirc_delt.ncml',...
    'dods://testbedapps.sura.org/thredds/dodsC/ugrid/TestCases/SELFE/selfe_delt.ncml'};

var='zeta';  % variable to plot

itimes=[30 8 280 15];  % which time step to read from each file

% There is nothing model specific in the loop below!
for i=1:length(uris)
    uri=uris{i};
    itime=itimes(i);
    tic
    %initialize NCUGRID dataset object
    nc=ncugrid(uri);
    nc.cells
    nc.nodes
    % read data at specified time step (all nodes)
    zeta = nc.data(var, [itime 1], [itime max(nc.size(var))], [1 1]);
    
    % read lon,lat & connectivity array!
    zgrid = nc.grid_interop(var, [itime 1], [itime max(nc.size(var))], [1 1]);
    
    lon = zgrid.lon;
    lat = zgrid.lat;
    ele = zgrid.connectivity;
    
    % plot the elevation
    figure(i)
    [m,n]=size(ele);if m==3,ele=ele.';end  % flip connectivity for trisurf
    trisurf(ele,lon,lat,zeta);shading interp;view(2);colorbar
    title(sprintf('%s (%s): %s',var,nc.attribute('units', var),datestr(zgrid.time)));
    set (gca, 'DataAspectRatio', [1 cos(mean(lat(find(isfinite(lat(:)))))*pi/180) 1] ); 
    toc
end
