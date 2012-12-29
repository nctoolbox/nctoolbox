function dat=nj_tseries(ncRef,var,loni,lati,varargin)
% NJ_TSERIES Interpolate time series at lon,lat points from CF model output
% Usage:  dat=nj_tseries(uri,var,lon,lat,varargin)
% Inputs:
%    ncRef = netcdf file name, dods url, or ncgeodataset object
%    var = variable name (string) (e.g. 'salt');
%    lon = vector of longitudes where time series are desired
%    lat = vector of latitudes  where time series are desired
% Optional additional parameter pairs
%    'method' =   'nearest' or 'linear'
%    'layer' =  number of layer to extract
%    'itime' =   indices of time to extract
%    'ele' =    connectivity array (required for non structured grids)
%
% Unstructured Grid Example:  surface layer salinity for 1st 100 time steps
% url='http://www.smast.umassd.edu:8080/thredds/dodsC/fvcom/hindcasts/30yr_gom3'
% dat=nj_tseries(url,'salinity',-70.5,42.5,'method','nearest',...
%     'layer',1,'itime',1:100,'ele',1);
%
% Structured Grid Example:  surface layer salinity for 1st 100 time steps
% url='http://geoport.whoi.edu/thredds/dodsC/coawst_2_2/fmrc/coawst_2_2_best.ncd'
% dat=nj_tseries(url,'salt',-70.5,42.5,'method','nearest',...
%     'layer',16,'itime',1:100);

% Outputs:
%      dat.vals = data values
%      dat.time = time values
%      dat.ii= i value
%      dat.jj= j value

% Rich Signell (rsignell@usgs.gov)

if (isa(ncRef, 'ncgeodataset')) %check for ncgeodataset Object
  nc = ncRef;
else
  % open CF-compliant NetCDF File as a Common Data Model (CDM) "Grid Dataset"
  nc = ncgeodataset(ncRef);
end
vart=nc.geovariable(var);
lon=nc.data(vart.getlonname);
lat=nc.data(vart.getlatname);

while ~isempty(varargin),
  switch lower(varargin{1}),
    case 'method',
      method=varargin{2};
    case 'layer',
      layer=varargin{2};
    case 'itime',
      itime=varargin{2};
    case 'ele',
      ele=varargin{2};
    otherwise,
      error(['Can''t understand property:' varargin{1}]);
  end
  varargin([1 2])=[];
end
if ~exist('method','var'),method='nearest';end
jdmat=nj_time(ncRef,var); % returns [] if no time variable exists
if ~exist('itime','var'),
  if ~isempty(jdmat)
    itime=1:length(jdmat);
  end
end

ntime=length(itime);
npoint=length(loni);
% initialize
u=zeros(ntime,npoint);
lono=zeros(npoint);
lato=zeros(npoint);
inear=zeros(npoint);
jnear=zeros(npoint);

switch method
  case 'nearest'  % find closest point
    for k=1:npoint
      if exist('ele','var');  % unstructured grid
        ind=nearxy(double(lon(:)),double(lat(:)),loni(k),lati(k));
        lono(k)=lon(ind);
        lato(k)=lat(ind);
        fprintf('extracting point %d:lon=%f,lat=%f\n',k,lono(k),lato(k));

        if ~exist('itime','var');
          u(:,k)=nc{var}(ind);  % no time variable, e.g. bathymetry
        elseif (exist('itime','var') && ~exist('layer','var'))
          u(:,k)=squeeze(nc{var}(itime,ind));
        elseif (exist('itime','var') &&  exist('layer','var'))
          u(:,k)=squeeze(nc{var}(itime,layer,ind));
          dat.layer=layer;
        else
          disp('invalid arguments');return
        end
        dat.ii(k)=ind;
      else   % structured grid
        if isvector(lon), [lon,lat]=meshgrid(lon,lat);end;
        ind=nearxy(double(lon(:)),double(lat(:)),loni(k),lati(k));
        [jj,ii]=ind2ij(double(lon),ind);
        lono(k)=lon(jj,ii);
        lato(k)=lat(jj,ii);
        fprintf('extracting point %d:lon=%f,lat=%f\n',k,lono(k),lato(k));
        if ~exist('itime','var');
          u(:,k)=nc{var}(jj,ii);  % no time variable, e.g. bathymetry
        elseif (exist('itime','var') && ~exist('layer','var'))
          u(:,k)=squeeze(nc{var}(itime,jj,ii));
        elseif (exist('itime','var') &&  exist('layer','var'))
          u(:,k)=squeeze(nc{var}(itime,layer,jj,ii));
          dat.layer=layer;
        else
          disp('invalid arguments');return
        end
        dat.ii(k)=ii;
        dat.jj(k)=jj;
      end
    end
  case 'linear'  % interpolate from surrounding points
    if exist('ele','var');  % unstructured grid
      error(' ''linear'' interpolation not working yet for unstructured grid. Use ''nearest'' ')
    else % structured grid
      lono=loni;
      lato=lati;
      if isvector(lon), [lon,lat]=meshgrid(lon,lat);end;
      [n,m]=size(lon);
      [ii,jj]=meshgrid(1:m,1:n);
      F1=TriScatteredInterp(double(lon(:)),double(lat(:)),ii(:));
      ivar=F1(loni(:),lati(:));
      F2=TriScatteredInterp(double(lon(:)),double(lat(:)),jj(:));
      jvar=F2(loni(:),lati(:));
      for k=1:length(loni);
        fprintf('interpolating at point %d:lon=%f,lat=%f\n',k,loni(k),lati(k));
        i0=floor(ivar(k));
        ifrac=ivar(k)-i0;
        j0=floor(jvar(k));
        jfrac=jvar(k)-j0;
        if ~exist('itime','var'); % 2D (no time dimension)
          uall=nc{var}(j0:j0+1,i0:i0+1);
          u1=uall(1,1);
          u2=uall(1,2);
          u3=uall(2,1);
          u4=uall(2,2);
        elseif (exist('itime','var') && ~exist('layer','var'))
          uall=squeeze(nc{var}(itime,j0:j0+1,i0:i0+1));
          u1=uall(:,1,1);
          u2=uall(:,1,2);
          u3=uall(:,2,1);
          u4=uall(:,2,2);
        elseif (exist('itime','var') &&  exist('layer','var'))
          uall=squeeze(nc{var}(itime,layer,j0:j0+1,i0:i0+1));
          u1=uall(:,1,1);
          u2=uall(:,1,2);
          u3=uall(:,2,1);
          u4=uall(:,2,2);
          dat.layer=layer;
        else               % 3D (2D with time dimension)
          disp('invalid arguments');return
        end
        dat.ii=ivar;
        dat.jj=jvar;
        ua=u1+(u2-u1)*ifrac;
        ub=u3+(u4-u3)*ifrac;
        u(:,k)=ua+jfrac*(ub-ua);
      end
    end
end
dat.vals=u;
dat.lon=lono;
dat.lat=lato;
if ~isempty(jdmat)
  dat.time=jdmat(itime);
end
