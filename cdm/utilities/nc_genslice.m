function [Tvar,Tdis,Tzed,Tlon,Tlat] = ...
  nc_genslice(file,varname,lonTrk,latTrk,timTrk,report)
% [Tvar,Tdis,Tzed,Tlon,Tlat] = roms_genslice(file,varname,lonTrk,latTrk,timTrk)
%
% Basically works like GSLICE. Get a "generic" vertical slice out of ROMS
% output file at coordinates specified by a vector of lon, lat and time
% (e.g. a ship transect or a glider track)
%
% Input:
%   file      netcdf or opendap url (aggregation across time)
%   varname   name of the variable we want along the track
%   lon, lat  coordinates along the track
%   timTrk    time values along the track
%                given as Matlab datenum
%                Must have the same dimension as lon/lat  OR  be a scalar in
%                which case it is assumed the track is for a single time
% Output:
%   Tvar      interpolated 2D data along the track
%   Tdis      2D range or distance along the track (in kilometers)
%   Tzed      2D depth along the track
%   Tlon      2D lon along the track
%   Tlat      2D lat along the track
%
%  gslice originally by John Evans
%  modified by Weifeng Gordon Zhang
%  This version by John Wilkin 2011-01-26
%  Alexander Crosby 2011 April - Update to using new toolbox

verbose = false;
ncks = false;
if nargin > 5
  if strcmp(report,'verbose')
    verbose = true;
  end
  if strcmp(report,'ncks')
    ncks = true;
    verbose = true;
  end
end

% instantiate dataset and variable - ACrosby
nc = ncgeodataset(file);
var = nc.geovariable(varname);

% Make all track coordinates inputs into vectors
lonTrk = lonTrk(:);
latTrk = latTrk(:);
timTrk = timTrk(:);
if length(timTrk)==1
  % assume we want the slice at this instantaneous time
  timTrk = timTrk*ones(size(lonTrk));
end
if ~isequal(size(lonTrk),size(latTrk),size(timTrk))
  error('The dimensions of lonTrk, latTrk and timTrk do not match')
end

% Determine lon/lat range that encompasses the track coordinates
ax = [min(lonTrk) max(lonTrk) min(latTrk) max(latTrk)];
st.lat = [ax(3) ax(4)];
st.lon = [ax(1) ax(2)];

% Determine the minimal model/data index range required to cover the track
% Iax,Jax will specify the subset of the full grid to read from netcdf
if verbose
  disp('Reading full domain lon/lat grid coordinates from ')
  disp([' ' file])
end

% [lonfull,latfull] = nj_lonlat(file,varname); ACrosby
g = var.grid_interop(:,:,:,:); % Assuming data is rank 4, time z X Y
lonfull = g.lon;
latfull = g.lat;
if isvector(lonfull)
  % catch case (not ROMS) where coordinates are 1-D vectors
  [lonfull,latfull] = meshgrid(lonfull,latfull);
end
[nJ,nI] = size(lonfull);
% [Jax,Iax] = lonlat2ij(lonfull,latfull,ax); ACrosby
[r1, r2, c1, c2] = geoij(var, st);
Jax = r1:r2;
Iax = c1:c2;

% This subset will be strictly inside the limits of the track, so
% expand Iax,Jax by 1 in each dimension (if possible) or some track
% points will fall outside the data subset
if Iax(1)~=0
  Iax = [Iax(1)-1 Iax];
end
if Jax(1)~=0
  Jax = [Jax(1)-1 Jax];
end
if Iax(end)~=nI
  Iax = [Iax Iax(end)+1];
end
if Jax(end)~=nJ
  Jax = [Jax Jax(end)+1];
end
disp('Spatial subset to be extracted from ')
disp([' ' file])
disp('  is for index limits ')
disp(['  J = ' int2str(Jax(1)) ',' int2str(Jax(end))])
disp(['  I = ' int2str(Iax(1)) ',' int2str(Iax(end))])

% Find the fractional time indices of the track times
if verbose
  disp('Reading full vector of all times from ')
  disp([' ' file])
end
% tim_mod = nj_time(file,varname); - ACrosby
tim_mod = g.time;
Ttrk = interp1(tim_mod,1:length(tim_mod),timTrk);

% ------------------------------------------------------------------------
if ncks
  % Give the ncks command that would extract the necessary subset
  % to a local file, then return
  % *** This only work for ROMS at present
  try
    modeltype = nc_attget(file,nc_global,'type');
    if strcmp(modeltype(1:4),'ROMS')
      coords = 's_rho,Cs_r,zeta,h,hc,';
      varlst = [' -v ' coords varname];
      deta = [' -d eta_rho,' int2str(Jax(1)) ',' int2str(Jax(end))];
      dxi = [' -d xi_rho,' int2str(Iax(1)) ',' int2str(Iax(end))];
      dt = [' -d ocean_time,' int2str(min(Ttrk)) ',' int2str(max(Ttrk)+1)];
      str = ['ncks -F' varlst dxi deta dt ' ' file ' '];
      disp(str)
      Tvar = str;
    else
      disp('You may need to modify dimension names for this model')
      Tvar = [];
    end
  catch
    disp('You may need to modify dimension names for this model')
    varlst = [' -v ' varname];
    deta = [' -d lat,' int2str(Jax(1)) ',' int2str(Jax(end))];
    dxi  = [' -d lon,' int2str(Iax(1)) ',' int2str(Iax(end))];
    dt   = [' -d time,' int2str(min(Ttrk)) ',' int2str(max(Ttrk)+1)];
    str  = ['ncks -F' varlst dxi deta dt ' ' file ' local.nc'];
    disp(str)
    Tvar = str;
  end
  return
end
% ------------------------------------------------------------------------

% Extract one time record just to get coordinates of the subset grid
if verbose
  disp('Reading subset lon/lat grid coordinates')
end
% [tmp,geo] = nj_grid_varget(file,varname,...
%   [1 1 Jax(1) Iax(1)],[1 Inf length(Jax) length(Iax)]);
geo = var.grid_interop(1, :, Jax(1):Jax(end), Iax(1):Iax(end));

% lon,lat of the subset grid
lon_mod = double(geo.lon);
lat_mod = double(geo.lat);
if isvector(lon_mod)
  [lon_mod,lat_mod] = meshgrid(lon_mod,lat_mod);
end

% I,J index into the subsetted grid - not the full grid
[I,J] = meshgrid(1:size(lon_mod,2),1:size(lon_mod,1));

% Load vertical coordinates into the structures that will hold the model
% output and times
if isvector(geo.z)
  % catch case (not ROMS) where coordinates are 1-D vectors
  geo.z = repmat(geo.z,[1 length(geo.lat) length(geo.lon)]);
end
geo_t1.z = -abs(geo.z);
geo_t2.z = -abs(geo.z);

% Get fractional I,J grid positions on the subset grid corresponding to
% the track locations
if verbose
  disp('Griddata is calculating i,j coordinates of the track ...')
end
Jtrk  = griddata(lon_mod,lat_mod,J,lonTrk,latTrk);
Itrk  = griddata(lon_mod,lat_mod,I,lonTrk,latTrk);
if verbose
  disp(' ... griddata done')
end

% Find track locations that are outside the model grid or available times
% Flag with valid==NaN. We will skip over these positions but leave NaNs in
% the output so that the space dimension of the input and output matches
valid = 1:length(lonTrk);
valid(isnan(Itrk)|isnan(Jtrk)|isnan(Ttrk)) = NaN;

if ~any(isfinite(valid))
  error('No track points are encompassed by the data. Check lon/lat/time')
end

% preallocate for speed
nz   = size(geo.z,1);
Tvar = NaN*ones([nz length(lonTrk)]);
Tzed = Tvar;
Tlon = Tvar;
Tlat = Tvar;
Tdis = Tvar;

% Alongtrack distance to be passed to output for plotting
dist = cumsum([0; sw_dist(latTrk,lonTrk,'km')]);

% Process the track in time chunks bracketed by model times
chunks = unique(floor(Ttrk));
chunks(isnan(chunks)) = [];

if verbose
  disp(['Processing this track requires ' int2str(length(chunks)) ...
    ' time intervals out of '])
  disp([' ' file])
  disp(' within time index limits ')
  disp([' T = ' int2str(min(Ttrk)) ':' int2str(max(Ttrk))])
end








% counter for alongtrack profiles
prof = 1;

% Open data object
% nc = mDataset(file); A Crosby : nc already open

for nchunk = 1:length(chunks)
  
  % index of time at beginning of interval
  tindex0 = chunks(nchunk);
  
  % find the chunk of track coordinates in this time interval
  Section = find(floor(Ttrk)==tindex0);
  lonS = lonTrk(Section);
  latS = latTrk(Section);
  timS = timTrk(Section);
  distS = dist(Section);
  JtrkS = Jtrk(Section);
  ItrkS = Itrk(Section);
  
  if verbose
    disp([' Doing interval ' int2str(nchunk) ' with ' ...
      int2str(length(Section)) ' coordinates in time range:'])
  end
  % Load time and data for this tindex and the next
  % We cannot economize by assuming that tindex increments by 1 for each
  % chunk because the track times might be at larger intervals than
  % the model output
%   geo_t1.time = nj_time(nc,varname,tindex0); ACrosby
  geo_t1.time = g.time(tindex0);
%   geo_t2.time = nj_time(nc,varname,tindex0+1); ACrosby
  geo_t2.time = g.time(tindex0+1);
%   data_t12 = nc{varname}(tindex0+(0:1),1:nz,Jax,Iax); ACrosby
  data_t12 = var.data(tindex0:(tindex0+1), 1:nz, Jax(1):Jax(end), Iax(1):Iax(end));
  data_t1 = squeeze(double(data_t12(1,:,:,:)));
  data_t2 = squeeze(double(data_t12(2,:,:,:)));
  
  if verbose
    disp(['  ' datestr(geo_t1.time) ' to ' datestr(geo_t2.time)])
  end
  
  for j = 1:length(Section)
    
    if verbose && rem(prof,10)==0
      disp(['   Profile ' int2str(prof) ' of ' int2str(length(lonTrk))])
    end
    
    if isfinite(valid(prof))
      % Interpolate to the track positions for both times using
      % simple linear weights in I,J fractional coordinates
      
      % space weights
      fJ = floor(JtrkS(j));
      fI = floor(ItrkS(j));
      wJ = JtrkS(j)-fJ;
      wI = ItrkS(j)-fI;
      Dwgt(1,1) = (1-wJ)*(1-wI);
      Dwgt(1,2) = (1-wJ)*wI;
      Dwgt(2,1) = wJ*(1-wI);
      Dwgt(2,2) = wJ*wI;
      
      % time weights
      Twgt_t2 = (timS(j)-geo_t1.time)/(geo_t2.time-geo_t1.time);
      Twgt_t1 = 1-Twgt_t2;
      
      % I have not allowed for partial weights in the event some of the
      % data are in the land mask
      
      % Get vertical vector of model values interpolated to
      % each track position at both times
      
      data = data_t1;
      Dslice_t1 = Dwgt(1,1)*data(:,fJ,fI) + Dwgt(1,2)*data(:,fJ,fI+1) + ...
        Dwgt(2,1)*data(:,fJ+1,fI) + Dwgt(2,2)*data(:,fJ+1,fI+1);
      data = geo_t1.z;
      Zslice_t1 = Dwgt(1,1)*data(:,fJ,fI) + Dwgt(1,2)*data(:,fJ,fI+1) + ...
        Dwgt(2,1)*data(:,fJ+1,fI) + Dwgt(2,2)*data(:,fJ+1,fI+1);
      
      data = data_t2;
      Dslice_t2 = Dwgt(1,1)*data(:,fJ,fI) + Dwgt(1,2)*data(:,fJ,fI+1) + ...
        Dwgt(2,1)*data(:,fJ+1,fI) + Dwgt(2,2)*data(:,fJ+1,fI+1);
      data = geo_t2.z;
      Zslice_t2 = Dwgt(1,1)*data(:,fJ,fI) + Dwgt(1,2)*data(:,fJ,fI+1) + ...
        Dwgt(2,1)*data(:,fJ+1,fI) + Dwgt(2,2)*data(:,fJ+1,fI+1);
      
      % Weight by fractional time to get final vertical vector
      Dslice = Twgt_t1*Dslice_t1 + Twgt_t2*Dslice_t2;
      Zslice = Twgt_t1*Zslice_t1 + Twgt_t2*Zslice_t2;
      
      % Enter this profile into the accumulated slice
      Tvar(:,prof) = Dslice;
      Tzed(:,prof) = Zslice;
      onez = ones([nz 1]);
      Tlat(:,prof) = latS(j)*onez;
      Tlon(:,prof) = lonS(j)*onez;
      Tdis(:,prof) = distS(j)*onez;
      
    else % position not valid
      
      % track coordinates are outside the data so enter NaNs for this point
      nanz = NaN*ones([nz 1]);
      Tvar(:,prof) = nanz;
      Tzed(:,prof) = nanz;
      Tlat(:,prof) = latS(j)*nanz;
      Tlon(:,prof) = lonS(j)*nanz;
      Tdis(:,prof) = distS(j)*nanz;
      
    end
    
    % update profile counter
    prof = prof + 1;
    
  end
  
end

% Close data object
% nc.close(); % Acrosby - not needed


