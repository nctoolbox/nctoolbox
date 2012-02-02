function [ U, g ] = cgrid_uv2rho(ncRef, uname, vname, hname, aname,...
    itime, klev, jj, ii)
% CGRID_UV2RHO Get 2D slice of velocity from a C-grid model (e.g. ROMS)
% Usage:
%    [U,g] = cgrid_uv2rho(ncRef, uname, vname, hvar, avar, itime, klev, [jj], [ii]);
% Inputs:
%   ncRef = OpenDAP Data URL, Local file, or existing 'ncgeodataset' object
%   uname - u variable to slice
%   vname - v variable to slice
%   hname - rho-centered variable to slice (e.g. 'h')
%   aname - rho-centered angle variable (e.g. 'angle')
%   itime = time step to extract data. Use -1 for last time step.
%   klev  = vertical level  Use -1 for last level.
%   [jj] =  jindex range (subset) [optional, default is entire dimension]
%   [ii] =  iindex range (subset) [optional, default is entire dimension]
%
% Outputs:
%   U = 2D array of complex velocity
%   g  = grid structure containing lon,lat,z,time (Matlab datenum)
%
% Example:
% url='http://geoport.whoi.edu/thredds/dodsC/examples/bora_feb.nc';
% hname='h'; uname='u'; vname='v'; aname='angle';
% nc=ncgeodataset(url);
% itime=3; % 3rd time step
% klev=-1; % last (top) layer
% Whole domain:
% [ U,g ] = cgrid_uv2rho(nc,uname,vname,hname,aname,itime,klev);
% Subset to jj,ii range:
% [ U,g ] = cgrid_uv2rho(nc,uname,vname,hname,aname,itime,klev,1:58,1:70);
% pcolorjw(g.lon,g.lat,abs(U));colorbar;dasp(44);
% arrows(g.lon(1:2:end,1:2:end),g.lat(1:2:end,1:2:end),...
%     U(1:2:end,1:2:end),0.08,'black');
%  title(datestr(g.time));dasp(44);

% Note: this routine assumes that the C grid has variables arranged:
%
%       ---------------------------------
%       rho | u | rho | u | rho | u | rho
%       ---------------------------------
%        v  |   | {v} |   | {v} |   |  v
%       ---------------------------------
%       rho |{u}|(rho)|{u}|(rho)|{u}| rho
%       ---------------------------------
%        v  |   | {v} |   | {v} |   |  v
%       ---------------------------------
%       rho | u | rho | u | rho | u | rho
%       ---------------------------------
%
% .. so that size(rho)=ny,nx;  size(u)=ny,nx-1, size(v)=ny-1,nx)
% Also assumes coordinate dimensions are arranged (t, z, y, x)

if (isa(ncRef, 'ncgeodataset')) %check for ncgeodataset object
    nc = ncRef;
else
    % open CF-compliant dataset as ncgeodataset object
    nc = ncgeodataset(ncRef);
end

uvar = nc.geovariable(uname);
vvar = nc.geovariable(vname);
hvar = nc.geovariable(hname);
avar = nc.geovariable(aname);

% depth from uvar
if nargin < 7
    usize = size(uvar);
    klev = 1:usize(2);
end

% get lon,lat size from hvar
if nargin < 8
    hsize = hvar.size;
    horder = hvar.getaxesorder;
    if horder.lat == horder.lon
        jj = 1:hsize(horder.lat);
        ii = 1:hsize(horder.lat  + 1);
    else
        %         ii = hsize(order.lat);
        %         jj = hsize(order.lon);
        error('order.lat != order.lon');
    end
else
    %pass
end
ujj = (jj(1) + 1):(jj(end) - 1);
uii = ii(1):ii(end) - 1;
vjj = jj(1):jj(end) - 1;
vii = (ii(1) + 1):(ii(end) - 1);

u = double(squeeze(uvar.data(itime, klev, ujj, uii))); % get u
v = double(squeeze(vvar.data(itime, klev, vjj, vii))); % get v
g = hvar.grid_interop(jj, ii); % get lon,lat from rho-point variable (like 'h' in ROMS)
ang = avar.data(jj, ii); % get angle
U = ones(size(g.lon)) * nan; % template for U at rho points
U(2:end-1, 2:end-1) = complex(av2(u.').', av2(v)); %average u,v to rho
U = U .* exp(sqrt(-1) * ang); % rotate

if nargout == 2
    tim = uvar.timewindowij( double(itime) );
    g.time = tim.time;
    g.klevel = klev;
    g.itime = itime;
end
