% COMP_UGRID_UV  
% Compare depth-averaged current from 2 unstructured grid
% models with UGRID conventions (http://bit.ly/cf_ugrid), allowing 
% comparison with no model specific code.  This code makes a movie
% if there is more than one time step specified

comp_titl='Comparison of 2D No-Wave runs for Hurricane Rita';

titl{1}='ADCIRC';
uris{1}='http://testbedapps.sura.org/threddsdev/dodsC/inundation/ADCIRC/ike/3Dvrwww'
%uris{1}='http://testbedapps.sura.org/threddsdev/dodsC/inundation/ADCIRC/rita/2Dvrwoww'
uname{1}='u-vel';
vname{1}='v-vel';
tname{1}='time';

titl{2}='FVCOM';
uris{2}='http://testbedapps.sura.org/threddsdev/dodsC/inundation/FVCOM/ike/3Dvrwww';
%uris{2}='http://testbedapps.sura.org/threddsdev/dodsC/inundation/FVCOM/rita/2Dvrwoww';
uname{2}='ua';
vname{2}='va';
tname{2}='time';

% for just one plot, define one time:
jdi=datenum([2008 9 13 03 0 0]); %Ike
%jdi=datenum([2005 9 24 06 0 0]); %Rita

% for a movie define more than one time
%jdi=datenum([2008 9 12 21 0 0]):3/24:datenum([2008 9 13 15 0 0]);
if length(jdi)>1
  movie_name='ike';
end

% scale factor for arrows
vfac=0.03;
% bounding box for figures
ax=[-95.2  -94.4   28.8   29.9]
% color range for figures
cax=[0 3];

if length(jdi)>1
  vidobj=VideoWriter(movie_name,'Uncompressed AVI');
  vidobj.FrameRate = 1;
  open(vidobj);
end

% There is nothing model specific below this point

%  loop through each model to read the grid information
for i=1:length(uris)
  tic
  %initialize dataset object
  nc{i}=ncgeodataset(uris{i});
  %get geovariable object
  uvar{i}=nc{i}.geovariable(uname{i});
  vvar{i}=nc{i}.geovariable(vname{i});
  jd{i}=nc{i}.time(tname{i});
  % get mesh variable name (inference array)
  gvar_name{i}=uvar{i}.attribute('mesh');
  if isempty(gvar_name{i}),
    gtype{i}='structured';
  else
    gtype{i}='unstructured';
  end
  switch gtype{i}
    case 'structured'
      g=uvar{i}.grid(itime,:,:);
      lon{i}=g.lon;
      lat{i}=g.lat;
    case 'unstructured'
      gridvar=nc{i}.variable(gvar_name{i});
      a=textscan(gridvar.attribute('node_coordinates'),'%s %s');
      lon{i}=nc{i}{a{1}}(:);
      lat{i}=nc{i}{a{2}}(:);
      location{i}=uvar{i}.attribute('location');
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
  end
end
% loop through time, reading data at each step from all models 
fig=figure(1);
clf;set(gcf,'renderer','zbuffer');set(gcf,'color','white');
for k=1:length(jdi);
  for i=1:length(uris)
    itime=date_index(jd{i},jdi(k));
    subplot(1,2,i);cla
    switch gtype{i}
      case 'structured'
        u=double(squeeze(uvar{i}.data(itime,:,:)));
        v=double(squeeze(vvar{i}.data(itime,:,:)));
        U=complex(u,v);
        pcolorjw(lon{i},lat{i},abs(U));shading interp;colorbar;...
          axis(ax);caxis(cax);
        arrows(lon{i},lat{i},U,0.03,'black');
      case 'unstructured'
        % read data at specified time step for all nodes
        u=double(squeeze(uvar{i}.data(itime,:)));
        v=double(squeeze(vvar{i}.data(itime,:)));
        U=complex(u,v);
        U(U==0)=nan;
        igood=isfinite(U);
        % handle case of 'face' location (e.g. FVCOM)
        if strmatch('face',location{i}),
          % calculate u,v lon,lat positions (faster than reading lonc,latc)
          lon_vel=1/3*(lon{i}(tri{i}(:,1))+lon{i}(tri{i}(:,2))+lon{i}(tri{i}(:,3)));
          lat_vel=1/3*(lat{i}(tri{i}(:,1))+lat{i}(tri{i}(:,2))+lat{i}(tri{i}(:,3)));
          lonp=lon{i}(tri{i});
          latp=lat{i}(tri{i});
          patch(lonp.',latp.',abs(U));view(2);shading flat % shading 'interp' doesn't work with patch
          arrows(lon_vel(igood),lat_vel(igood),U(igood),vfac,'black');
        else
          trisurf(tri{i},lon{i},lat{i},0*lon{i}-5,abs(U));view(2);shading interp;
          arrows(lon{i}(igood),lat{i}(igood),U(igood),vfac,'black');
        end
        colorbar;axis(ax);caxis(cax);...
          fac=cos(mean(lat{1}(:))*pi/180);...
          title(sprintf('%s %s (%s): %sZ',titl{i},uname{i},...
          uvar{i}.attribute('units'),datestr(jd{i}(itime),'yyyy-mm-dd HH:MM')));...
          set (gca, 'DataAspectRatio', [1 fac 1] );...
          set (gca, 'tickdir','out');
        toc
        
    end
  end
  if(length(jdi)>1)
    F=getframe(fig);
    writeVideo(vidobj,F);
  end
end
if(length(jdi)>1)
  close(vidobj)
end