% DEMO2A - Subsetting data with Matlab style indexing

echo('on')
% STARTING DEMO2A ----------------------------------------------------------
% An example of subsetting data using ncdataset

%% ---- Open the dataset
ds = ncgeodataset('http://geoport.whoi.edu/thredds/dodsC/examples/OS_M1_20081008_TS.nc')

% You can view the variables available to you
ds.variables

%% Plot all the data
figure;
plot(ds.time('TIME'), ds{'TEMP'}(1:max(ds.size('TEMP')), 1, 1, 1))
hold('on')

%% ---- Lets fetch a subset of time in Matlab's native format
startIdx = 100;
endIdx = max(ds.size('TIME'));
stride = 10;
t = ds.data('TIME', startIdx, endIdx, stride);
t = ds.time('TIME', t); % Convert time data to matlab format. See help ncdataset.time

%% ---- Now lets get a subset of the temperature data.
% NOTE: The shape of the variables size is important for subsetting
ds.size('TEMP')
% Use variable, start, end, stride to subset
% Same as: temp = ds.data('TEMP', [startIdx 1 1 1], [endIdx 1 1 1], [stride 1 1 1]);
temp = ds{'TEMP'}(startIdx:stride:endIdx, 1, 1, 1);


%% ---- Add Subsetted Data to Plot 
plot(t, temp, 'r.');...
datetick('x', 2);grid;...
legend('All Data', 'Decimated Data');...
title('Surface Temperature at M1 Mooring in Monterey Bay');...
ylabel('Temperature [^oC]');
hold('off');
shg

echo('off') % ENDING DEMO2A ------------------------------------------------