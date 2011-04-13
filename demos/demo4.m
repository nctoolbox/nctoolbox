% DEMO4 - Playing with NetCDF-JAVA
echo('on')
% Starting DEMO4 ----------------------------------------------------------
% Example of data access using the NetCDF-Java API
ds = ncdataset('http://elvis.shore.mbari.org/thredds/dodsC/all/OS_M1_R/OS_MBARI-M1_20091020_R_TS.nc');
nc = ds.netcdf;

%% ---- Query NetCDF using Java API
nc.getConventionUsed()
nc.getDetailInfo()
nc.hasUnlimitedDimension()
nc.getUseNaNs() 

%% ---- Set properties using Java API
% The following line turns off the auto conversion of missing_values to NaN
% NOTE: This switches ALL netcdf files not just the one you are working on.
% nc.setUseNaNs(0)

%% ---- Dump NcML using Java API
out = java.io.BufferedWriter(java.io.FileWriter('trash_me.xml'));
nc.writeNcML(out, '')
out.close()
echo('off') % Ending DEMO4 ------------------------------------------------