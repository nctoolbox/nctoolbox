% DEMO3C - Chaining files using ncml and finding attributes using 'value4key'

echo('on')
o2_ = 'Oxygen';
time_ = 'esecs';

%% ---- Create the full pathname for an ncml file

scriptname = which('demo3c');  %% find the absolute path of a local file
fName = strrep(scriptname,'.m','.ncml') ; % create a path to the associated ncml file

if ~exist(fName,'file')
  error(strcat(fName,' does not exist.'));
end

%% ---- Access the OpenDAP datasets through an ncml file

nc = ncdataset(fName); % read the local ncml file
t = [nc.time(time_)];
o2 = [nc.data(o2_)];

%% ---- Find the units of Oxygen and label the display
attr = nc.attributes(o2_);
units = value4key(attr, 'units'); % Retrieve the units value
name = value4key(attr, 'long_name'); % Retrieve the long_name value

%% ---- Plot the data
figure;
plot(t, o2);...
datetick('x');...
ylabel([name ' [' units ']']);...
title({'M1 Mooring in Monterey Bay',nc.location},'interpreter','none');...
grid;...
shg
echo('off')
