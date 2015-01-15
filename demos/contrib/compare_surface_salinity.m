% COMPARE_SURFACE_SALINITY
% compare surface salinity from various models:
% 1. query metadata database for datasets in specified space/time range
% 2. get data using OPeNDAP endpoints discovered in metadata
% 3. Use NCTOOOLBOX to allow interoperable geospatial subsetting of data
any_text='sea_water_salinity';
csw_endpoint = 'http://www.ngdc.noaa.gov/geoportal/csw';
scheme='urn:x-esri:specification:ServiceType:odp:url';
bbox=[-75.0, -71.0, 39.0, 41.5];   % lon_min lon_max lat_min lat_max
%bbox=[ -76.5220,  -75.7105,   36.8248,   37.7850];
s.bbox=bbox;
s.start=now-3;
s.stop=now+3;
s.any_text=any_text;
s.endpoint=csw_endpoint;
cax=[32.0 35.0];
[records,paramString] = csw_search(s);
scheme='urn:x-esri:specification:ServiceType:odp:url';
datasets = csw_scheme(records,scheme);
%%
k=0;  % figure number
for i=1:length(datasets)
    % choose only datasets that don't have 'Averages' in the title.
    if isempty(strfind(datasets{i}.title,'FVCOM')),
        disp(datasets{i}.title)
        % open as geodataset using NCTOOLBOX
        nc=ncgeodataset(datasets{i}.url);
        % find the variable that has the specified standard_name
        var = find_std_names(nc,any_text);
        ncgvar=nc.geovariable(var);
        % subset using geo coordinates, not indices!
        s.time=datenum(datestr(s.start));
        s.lon=bbox([1 2]);
        s.lat=bbox([3 4]);
        % here the only hack, because in roms native, level[1] is bottom
        % we should add a feature to geosubset to specify 'top' or 'bottom'
        %       if strfind(datasets{i}.title,'ROMS'),
        s.z_index='surface';
        %       else
        %           s.z_index=1;
        %       end
        % adding try/catch here to catch cases where geosubset doesn't find
        % any data within the specified region.  This can happen when the
        % metadata bounding box is bigger than the region of valid data.
        try
            sub = ncgvar.geosubset(s);
            
            % just plotting stuff below here
            k=k+1
            figure(k);
            pcolorjw(sub.grid.lon,sub.grid.lat,double(squeeze(sub.data)));
            set(gca,'tickdir','out');
            set(gcf,'color','white');
            set(gca,'xgrid','on','ygrid','on','layer','top');
            set(gcf,'color',[0.85 0.85 0.85])
            caxis(cax);
            colorbar;
            title([datasets{i}.title ':' any_text ':' datestr(s.time)],'interpreter','none')
        catch
        end
    end
end