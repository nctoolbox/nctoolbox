function [dataout] = ncstationlook(varargin)

inputs = arg2hash(varargin);
obs = value4key(inputs, 'observation');
obsvar = value4key(inputs, 'obs_variable');
models = value4key(inputs, 'model');
modvar = value4key(inputs, 'mod_variable');
modvar2 = value4key(inputs, 'mod_variable2');
date = value4key(inputs, 'date');

if ~iscell(obs)
    obs = {obs};
end
if ~iscell(models)
    models = {models};
end

%
% figure
% usa = shaperead('usastatelo');
% mapshow(usa);
%  hold on;
tic
for i = 1:length(models)
    nc_mod = ncgeodataset(models{i});
    try
        modvariable = nc_mod.geovariable(modvar);
    catch
        modvariable = nc_mod.geovariable(modvar2);
    end
    
    
    for j = 1:length(obs)

        nc_obs = ncgeodataset(obs{j});
        
        obsvariable = nc_obs.geovariable(obsvar);
        
        
        obstime = obsvariable.timewindow(date);
        modtime = modvariable.timewindow(obstime.time);
        
        obsdata = obsvariable.data(obstime.index, :);
        moddata = modvariable.data(modtime.index, :);
        
        obsgrid = obsvariable.grid_interop(obstime.index, :);
        modgrid = modvariable.grid_interop(modtime.index, :);
        
        modstruct.var_struct.grid = modgrid;
        modstruct.var_struct.data = moddata;
        
        [modprofile lat lon] = profilefrom4d(modstruct, [obsgrid.lon obsgrid.lat]);
        obsprofile = obsdata(~isnan(obsdata));
        obsdepths = obsgrid.z(~isnan(obsdata));
        
        modprofile.depths = ncunits(modprofile.depths,nc_mod.attribute('Grid_Depth','units'), 'm');
        obsdepths = ncunits(obsdepths,'ft', 'm');
        
        h = figure;
       
        b = subplot(121);
        plot(obsprofile, obsdepths);
        hold on
        scatter(obsprofile, obsdepths);
        title({[nc_obs.attribute('title')];['at ', num2str(obsgrid.lon), ', ', num2str(obsgrid.lat)];['on ', num2str(datestr(obstime.time))]})
        xlabel([obsvariable.name,' in ', obsvariable.attribute('units') ])
        ylabel(['Depth in m'])
        
        v = subplot(122);
        plot(modprofile.profile, modprofile.depths);
        hold on
        scatter(modprofile.profile, modprofile.depths);
        title({[nc_mod.attribute('title')];['at ', num2str(lon),', ', num2str(lat)];['on ', num2str(datestr(modtime.time))]})
        xlabel([modvariable.name,' in ', modvariable.attribute('units') ])
        ylabel(['Depth in m'])
        
        stand = std([modprofile.profile, obsprofile]);
        xlim([min([modprofile.profile, obsprofile])-stand max([modprofile.profile, obsprofile])+stand]);
        ylim([(min([modprofile.depths, obsdepths])-5) 0])
        linkaxes([b v]);
        toc
    end
end