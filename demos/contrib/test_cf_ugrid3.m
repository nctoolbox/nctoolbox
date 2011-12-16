% TEST_CF_UGRID3
% Read/plot CF-UGRID Convention data from ADCIRC and SELFE
% next up: FVCOM and ELCIRC

titl{1}='ADCIRC';
uris{1}='http://testbedapps.sura.org/threddsdev/dodsC/inundation/ADCIRC/ike/3Dvrwww'
vars{1}='zeta';
times{1}=[2008 9 13 06 0 0];

titl{2}='SELFE';
uris{2}='http://testbedapps.sura.org/threddsdev/dodsC/inundation/selfe/ike/3Dvrwww';
vars{2}='elev';
times{2}=[2008 9 13 06 0 0];

% bounding box for figures
ax=[-95.4519  -87.3856   27.3950   31.3868]
% color range for figures 
cax=[0 5];

% There is nothing model specific in the loop below!
for i=1:length(uris)
    tic
    %initialize dataset object
    nc=ncgeodataset(uris{i});
    %get geovariable object
    zvar=nc.geovariable(vars{i});
    % Find the coordinate variables
    lon=nc{zvar.getlonname}(:);
    lat=nc{zvar.getlatname}(:);
    tdat=zvar.timewindowij(times{i});
    itime=tdat.index;
    % read data at specified time step for all nodes
    zeta=zvar.data(itime,:);

    % get grid variable name (inference array)
    gvar=zvar.attribute('mesh');

    % get data from grid variable (connectivity array)
    grid=nc{gvar}(:);
    [m,n]=size(grid);
    % check/fix orientation of connectivity array
    if m==3,
        grid=grid.';
    elseif n~=3
        disp('Error:Currently handling triangles only');return
    end
    figure(i)
    trisurf(grid,lon,lat,zeta);shading interp;view(2);colorbar;...
    axis(ax);caxis(cax);...
    title(sprintf('%s %s (%s): %sZ',titl{i},vars{i},...
      zvar.attribute('units'),datestr(tdat.time)));...
    % set aspect ratio for lon/lat plot based on mean latitude
    set (gca, 'DataAspectRatio', [1 cos(mean(lat(:))*pi/180) 1] );

    toc
end
