% Script to access surface currents from Fukushima ocean forecast
% via OPeNDAP using NCTOOLBOX 
% http://code.google.com/p/nctoolbox/

%  This script uses latest version from Mecurial, which contain
%  a new matlab-based indexing not yet in the public release:
%  hg clone https://nctoolbox.googlecode.com/hg/ nctoolbox 
%  Rich Signell (rsignell@usgs.gov)

url='http://edac-dap3.northerngulfinstitute.org/thredds/dodsC/ncom_fukushima_agg/Fukushima_best.ncd';
nc=cfdataset(url);
dn=nc.time('time');
u=nc.variable('water_u');
v=nc.variable('water_v');
ksurface=1;
%isub=1;
isub=4;  %optionally subsample lon/lat array for speed/testing
sz=size(u);
grd=u.grid(1,1,1:isub:sz(3),1:isub:sz(4)); %nb: "1:end" does not work yet
for itime=1:length(dn);
    udata=u.data(itime,ksurface,1:isub:sz(3),1:isub:sz(4));
    vdata=v.data(itime,ksurface,1:isub:sz(3),1:isub:sz(4));
    disp(sprintf('Step %d/%d',itime,length(dn)))
    w=double(squeeze(complex(udata,vdata)));
    pcolor(grd.lon,grd.lat,abs(w));shading flat;colorbar
    title(sprintf('NCOM Fukushima forecast surface currents (m/s): %s',...
        datestr(dn(itime))));
    caxis([0 2]);set(gca,'tickdir','out');shg
end