% DEMO5 - Accessing local and webserver data (non-OpenDAP)

echo('on')
% Starting DEMO5 ----------------------------------------------------------
% Accessing both local and remote NetCDF files.
%% ---- Open local and remote datasets
ds = {  ...
    ncdataset('m1_pco2.nc'), ... % ---- Open local NetCDF
    ncdataset('http://www3.mbari.org/staff/brian/pub/m1_pco2.nc') ... % ---- Open NetCDF on webserver
};

%% ---- Plot them together
p = 'pco2';
t = 'time_d';
colors = {'r.', 'g.'};
for i = 1:length(ds)
    plot(ds{i}.time(t), ds{i}.data(p), colors{i});
    hold('on')
end
hold('off')
datetick('x')
set(gca, 'XLim', datenum(['2005-01-01'; '2008-01-01']))
grid('on')
legend('local NetCDF', 'NetCDF on web server')
text(.05,.9,{ds{1}.location, ds{2}.location},'Units','normalized','interpreter','none')
a = ds{1}.attributes(p)
ylabel([value4key(a, 'long_name') ' [' value4key(a, 'units') ']'])
title('pCO_2 at M1 Mooring in Monterey Bay')
echo('off') % Ending DEMO5 ------------------------------------------------


