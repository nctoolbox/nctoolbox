% DEMO10 - Example of working with Java Structures using Netcdf-Java API. 
% This is just a reference until we decide how to best wrap this in a Matlab API.
% Also, this code is probably naive, I don't usually work with Structures. The 
% samples represent CTD bottle data.

% Brian Schlining
% 2010-02-02

echo('on')
% Starting DEMO10 ---------------------------------------------------------
% NCTOOLBOX doesn't yet support structures directly. But you can work with using
% the netcdf-Java API in Matlab


ds = ncdataset('http://elvis.shore.mbari.org/dods/bog/BCTD');
% URL inaccessible with 404 2014-04-03

%% netcdf.getVariables() returns a Java ArrayList we use the java method get(0)
% to return the first (and in this case only) variable/structure. The result from
% 'get' is also a Java ArrayList
bctd = ds.netcdf.getVariables().get(0)

%% Use NetCDF-Java copyto*JavaArrary()
% Read in the data of interest. I'm just reading in all the data for each variable
% of interest, then subsetting it in memory.
depth = bctd.findVariable('DEPTH').read().copyToNDJavaArray();
deep = find(depth > 4300);
lat = bctd.findVariable('DEC_LAT').read().copyToNDJavaArray();
lon = bctd.findVariable('DEC_LONG').read().copyToNDJavaArray();
siteId = bctd.findVariable('SITE_ID').read().copyTo1DJavaArray();
tmp = bctd.findVariable('TMP').read().copyTo1DJavaArray();
o2 = bctd.findVariable('O2').read().copyToNDJavaArray();

echo('off') % Ending DEMO10 -----------------------------------------------
% Dump the results to screen
fprintf(1, 'DEPTH, DEC_LAT, DEC_LONG, SITE_ID, TMP, O2\n');
for i = 1:length(deep)
    j = deep(i);
    fprintf(1, '%4i, %7.3f, %7.3f, %s, %6.4f, %5.3f\n', depth(j), lat(j), lon(j), ...
            siteId(j), tmp(j), o2(j));
end

