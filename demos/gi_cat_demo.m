% Plot all time series variables that have a certain
% 'standard_name' in a specified bounding box and time period

open_url='http://geoport.whoi.edu/gi-cat/services/opensearch';
var='sea_water_temperature'; %free text search string
time_var='time'; %time coordinate variable
bbox =[-180 180 0 90];  % [lon_min lon_max lat_min lat_max]
start=[1990 1 1 0 0 0]; %start
stop =[2000 1 1 0 0 0];  %stop

% opensearch query
q.endpoint=open_url;
q.bbox=sprintf('%d,%d,%d,%d',bbox([1 3 2 4]));
q.time_start=datestr(start,'yyyy-mm-ddTHH:MM:SSZ')% convert to ISO
q.time_end=datestr(stop,'yyyy-mm-ddTHH:MM:SSZ')% convert to ISO
q.string_text=var;
[links,params]=opensearch(q);   % make the query

dap=links2dap(links); % find only the OPeNDAP links

for i=1:length(dap);
    figure(i);
    nc=cfdataset(dap{i});
    vars=nc.variables;
    for j=1:length(vars); %loop through variables to find standard_names
        std_name=value4key(nc.attributes(vars{j}),'standard_name');
        if strcmp(std_name,var),
            vart=nc.variable(vars{j});
            jd=nc.time(time_var);
            ii=date_index(jd,start,stop);  % find indices of dates between start/stop
            jd=jd(ii);
            t=vart.data(ii);  %extact these indices from dataset
            plot(jd,t);
            ylabel(sprintf('%s [%s]',var,value4key(vart.attributes,'units')),...
                'interpreter','none');
            datetick
            grid;
            title(value4key(nc.attributes,'title'))
        end
    end
end