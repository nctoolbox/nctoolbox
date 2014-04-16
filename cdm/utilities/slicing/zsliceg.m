function varz = zsliceg(data, z, depth)
% ZSLICEG: horizontal slice 3D data at fixed vertical level
%         using generic input structures (e.g returned by the NJ Toolbox)
% Usage: varz=zsliceg(data,z,depth);
% Input: data = 3D array of data
%        z    = 1D or 3D array of z values (if 1d, will be expanded to 3d)
%        depth = desired z level for horizontal slice
% Output: varz = 2D array of data
% Example:
% uri='http://geoport.whoi.edu/thredds/dodsC/examples/bora_feb.nc'
% [t,g]=nj_tslice(uri,'temp',1);% grab 3d field of 'temp' at time step 1
% tz = zsliceg(t,g.z,-5);  % return temperature slice at 5 m depth
% pcolorjw(g.lon,g.lat,double(tz));  % plot temp at 5 m depth
%
% NCTOOLBOX (https://github.com/nctoolbox/nctoolbox)

[NZ,NY,NX]=size(data);
[nz,ny,nx]=size(z);
if(nx==1 & ny==1),
  z3d(:,1,1)=z;
  z=repmat(z3d,[1 NY NX]);
end
zdiff=z(1,:,:)-z(end,:,:);
if min(zdiff(:)) < 0;  % if z(end) is top, then flipud so that z(1)=top
  z = flipdim(z,1);
  data = flipdim(data,1);
end
iland=find(~isfinite(data(1,:)));
z = reshape(z,[NZ NX*NY]);
data = reshape(data,[NZ NX*NY]);
ind=find(isfinite(z(1,:)));


% Add a layer of z values above surface with Inf
% Add a layer of z values below bottom with -Inf
z = [Inf*ones(1,NX*NY); z; -Inf*ones(1,NX*NY)]; 

% Add a layer of data values above surface cloning the surface values
% Add a layer of data values below bottom with Nan
data = [data(1,:); data; NaN*ones(1,NX*NY)];

% Find the indices of data values that have just greater depth than
% depth

zg_ind = find(diff(z<depth)~=0);
zg_ind = zg_ind + [0:1:length(zg_ind)-1].';

data_greater_z = data(zg_ind);
depth_greater_z = z(zg_ind);

% Find the indices of the data values that have just lesser depth
% than depth
zl_ind = find(diff(z>depth)~=0);
zl_ind = zl_ind + [1:1:length(zg_ind)].';
data_lesser_z = data(zl_ind);
depth_lesser_z = z(zl_ind);

% Linearly interpolate between the data values.
alpha = (depth-depth_greater_z)./(depth_lesser_z-depth_greater_z);
data_at_depth = (data_lesser_z.*alpha)+(data_greater_z.*(1-alpha));
varz = reshape(data_at_depth,[NY NX]);
varz(iland)=nan;  % mask out land

end
