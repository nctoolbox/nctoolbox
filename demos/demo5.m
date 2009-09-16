% DEMO5 - Accessing local and webserver data (non-OpenDAP)

echo('on')
%% ---- Open datasets
ds = {  ...
    ncdataset('m1_pco2.nc'), ... % ---- Open local NetCDF                                 
    ncdataset('http://www.mbari.org/staff/brian/pub/m1_pco2.nc') ... % ---- Open NetCDF on webserver
}; 

%% ---- Plot them together
p = 'pco2';
t = 'time_d';
colors = {'r.', 'g.'};
for i = 1:length(ds)
    plot(ds{i}.time(t), ds{i}.data(p), colors{i});
    hold('on')
end
datetick('x')
grid('on')
legend('local NetCDF', 'NetCDF on web server')
echo('off')


