% TEST_CF_UGRID2
% Compare Unstructured & Structured Grid models in N. Gulf of Mexico
clear
figure(1)

i=1;
titl{i}='NGOFS'; %FVCOM
grid_type{i}='unstructured_grid';
uris{i}='http://geoport-dev.whoi.edu/thredds/dodsC/ngofs/forecast';
vars{i}='zeta';
times{i}=[2012 1 6 12 0 0];

% i=i+1;
% titl{i}='SELFE';
% uris{i}='http://testbedapps-dev.sura.org/thredds/dodsC/inundation/selfe/ike/3Dvrwww';
% vars{i}='elev';
% times{i}=[2008 9 13 06 0 0];
% grid_type{1}='unstructured_grid';

i=i+1;
titl{i}='AMSEAS'; %NCOM
grid_type{i}='structured_grid';
uris{i}='http://edac-dap3.northerngulfinstitute.org/thredds/dodsC/ncom_amseas_agg/AmSeas_Aggregation_best.ncd';
vars{i}='surf_el';
times{i}=[2012 1 6 12 0 0];



% bounding box for figures
ax=[-98  -85.0   25.5   31.0];
% color range for figures
cax=[-.6 .0];
% END OF USER SPECIFIED PARAMETERS

% There is nothing model specific in the loop below!
for i = 1:length(uris)
  tic
  %initialize dataset object
  nc{i}=ncgeodataset(uris{i});
  %get geovariable object
  zvar{i}=nc{i}.geovariable(vars{i});
  % Find the coordinate variables
  lon{i}=nc{i}{zvar{i}.getlonname}(:);
  lon{i}(lon{i}>180)=lon{i}(lon{i}>180)-360; % convert lon to [-180,180] range
  lat{i}=nc{i}{zvar{i}.getlatname}(:);
  tdat=zvar{i}.timewindowij(times{i});
  itime=tdat.index;
  subplot(length(uris),1,i)
  jd{i}=tdat.time;
  switch grid_type{i}
    case 'unstructured_grid',
      % read data at specified time step for all nodes
      zeta{i}=zvar{i}.data(itime,:);
      
      % get mesh variable name (inference array)
      gvar_name=zvar{i}.attribute('mesh');
      gridvar=nc{i}.geovariable(gvar_name);
      % find connectivity array variable from mesh variable
      trivar_name=gridvar.attribute('face_node_connectivity');
      
      % get connnectivity array data
      tri{i}=nc{i}{trivar_name}(:);
      [m,n]=size(tri{i});
      % check/fix orientation of connectivity array
      if m==3,
        tri{i}=tri{i}.';
      elseif n~=3
        disp('Error:Currently handling triangles only');return
      end
      trisurf(tri{i},lon{i},lat{i},zeta{i});shading interp;view(2);...
        colorbar;axis(ax);caxis(cax);
    case 'structured_grid',
      s.lon=[ax(1) ax(2)];
      s.lat=[ax(3) ax(4)];
      s.t_index=itime;
      sub=zvar{i}.geosubset(s);
      lon{i}=sub.grid.lon;lat{i}=sub.grid.lat;
      lon{i}(lon{i}>180)=lon{i}(lon{i}>180)-360;
      zeta{i}=double(squeeze(sub.data));
      jd{i}=sub.grid.time;
      pcolorjw(lon{i},lat{i},zeta{i});
      colorbar;axis(ax);caxis(cax)
  end
  units{i}=zvar{i}.attribute('units');
  title(sprintf('%s %s (%s): %sZ',titl{i},vars{i},units{i},...
    datestr(jd{i})),'interpreter','none');...
    % set aspect ratio for lon/lat plot based on mean latitude
  set (gca, 'DataAspectRatio', [1 cos(mean(lat{1}(:))*pi/180) 1] );
  
  toc
end
set(gcf,'color','white');

figure(2)
%comp_ngulf
%cax=[0 .3]
 [x,y]=meshgrid(lon{2},lat{2});
 F=TriScatteredInterp(x(:),y(:),double(zeta{2}(:)));
 zetai=F(double(lon{1}(:)),double(lat{1}(:)));
 subplot(211)
 trisurf(tri{1},lon{1},lat{1},zeta{1});shading interp;view(2);...
        colorbar;axis(ax);caxis(cax);
        title(sprintf('%s %s (%s): %sZ',titl{1},vars{1},units{1},...
    datestr(jd{1})),'interpreter','none');...
    % set aspect ratio for lon/lat plot based on mean latitude
  set (gca, 'DataAspectRatio', [1 cos(mean(lat{1}(:))*pi/180) 1] );
  
 subplot(212)
 zoff=mean(zeta{1}(:))-mean(zetai(isfinite(zetai))); %difference in mean
 trisurf(tri{1},lon{1},lat{1},zetai+zoff);shading interp;view(2);...
        colorbar;axis(ax);caxis(cax);
        title(sprintf('%s %s (%s): %sZ',titl{2},vars{2},units{2},...
    datestr(jd{1})),'interpreter','none');...
    % set aspect ratio for lon/lat plot based on mean latitude
  set (gca, 'DataAspectRatio', [1 cos(mean(lat{1}(:))*pi/180) 1] );
  
  set(gcf,'color','white');

disp('Click left button to draw, right button to quit')
ibutton=1
while ibutton==1
   figure(2);
   [xx,yy,ibutton]=ginput(1);
  if(ibutton==1),
   line(xx,yy, 'color','black','Marker','o');
   dat1=nj_tseries(nc{1},'zeta',xx,yy,'method','nearest','ele',1);
   tdat=zvar{2}.timewindowij(dat1.time(1),dat1.time(end));
   dat2=nj_tseries(nc{2},'surf_el',xx+360,yy,'method','nearest','itime',tdat.index);
   figure(3);
   plot(dat1.time,dat1.vals,dat2.time,dat2.vals+zoff);legend('NGOFS','AMSEAS')
   datetick;grid
  end
  
end
