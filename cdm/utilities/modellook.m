% MODELLOOK - Infrastructure for quick model visualization.
%
% Useful to look at 2d, 3d, and 4d variables from gridded model netcdf files.
%
% Usage : >> modellook(daplink, datevec)
%             >> modellook(daplink, datenum)
%
% NCTOOLBOX (http://code.google.com/p/nctoolbox)
%
function [] = modellook(dap, date)


nc = ncgeodataset(dap);
vars = nc.variables;
for i = 1:length(vars)
    sizes{i} = [vars{i}, mat2str(nc.size(vars{i}))];
   
end

[select, ok] = listdlg('ListString', sizes,'SelectionMode','single', 'Name', 'ModelLook: Quick nc vis',...
    'PromptString', 'Please Select Variable:' );
var = nc.geovariable(vars{select});
try
    tind = var.timewindow(date);
    var_struct.data = var.data(tind.index, :, :, :);
    var_struct.grid = var.grid_interop(tind.index, :, :, :);
    out.ind = tind.index;
catch
    switch length(var.size)
        case 3
    var_struct.data = var.data(1, :, :);
    var_struct.grid = var.grid_interop(1, :, :);
        case 2
            var_struct.data = var.data(:, :);
    var_struct.grid = var.grid_interop(:, :);
    end
end
a = @press;
figure;
out.var = var;
out.var_struct = var_struct;


set(gcf, 'UserData', out);


if length(var.size) > 3
    h = pcolor(squeeze(double(var_struct.grid.lon)), squeeze(double(var_struct.grid.lat)),...
        squeeze(double(var_struct.data(1, 1, :,:)))); % remember to change back
    title({'ModelLook:';'Ctl-Click to create a path for vertical section and Shift-Click to plot.';...
        'Double-Click to create a profile plot at model location.';...
        ' '; [var.name, ' in ', var.attribute('units'), ' on ', datestr(tind.time)]}, 'interpreter', 'none');
elseif length(var.size) > 2
    if isfield(var_struct.grid, 'z')
        h = pcolor(squeeze(double(var_struct.grid.lon)), squeeze(double(var_struct.grid.lat)),...
            squeeze(double(var_struct.data(1, :,:)))); 
        title({'ModelLook:';...
            [var.name, ' in ', var.attribute('units')]}, 'interpreter', 'none');
    elseif isfield(var_struct.grid, 'time')
        h = pcolor(squeeze(double(var_struct.grid.lon)), squeeze(double(var_struct.grid.lat)),...
            squeeze(double(var_struct.data(1, :,:)))); 
        title({'ModelLook:';...
            'Double-Click to create a timeseries plot at model location.';...
            ' '; [var.name, ' in ', var.attribute('units')]}, 'interpreter', 'none');
    else
    end
elseif length(var.size) > 1
    h = pcolor(squeeze(double(var_struct.grid.lon)), squeeze(double(var_struct.grid.lat)),...
        squeeze(double(var_struct.data(:,:))));
    title({'ModelLook:';...
        ' '; [var.name, ' in ', var.attribute('units')]}, 'interpreter', 'none');
end
set(h, 'ButtonDownFcn', a)
shading interp
set(h, 'EdgeColor', 'k', 'EdgeAlpha', .3);

colorbar
xlabel(nc.attribute('title'))


    function [out temp tspoint] = press(src, evnt)
        sel_typ = get(gcbf,'SelectionType');
        switch sel_typ
            case 'alt'
                out = get(gcf, 'UserData');
                temp = get(gca,'CurrentPoint');
                if ~isfield(out, 'path')
                    out.path = [];
                end
                out.path(end+1,1:2) = temp(1,1:2);
                set(gcf, 'UserData', out);
            case 'open'
                temp = get(gcf, 'UserData');
                temp.points = get(gca,'CurrentPoint');
                temp.tspoints = temp.points(1,1:2);
                
                temp.size = temp.var.size;
                    temp.data = [];
                    temp.depths = [];
                    temp.profile = [];
                    switch length(temp.size)
                        case 4
                            if length(temp.var_struct.grid.lon) ~= length(temp.var_struct.grid.lat);
                                [temp.var_struct.grid.lon temp.var_struct.grid.lat] = meshgrid(temp.var_struct.grid.lon, temp.var_struct.grid.lat);
                            else
                            end
                            if length(size(temp.var_struct.grid.z)) > 2
                                
                                temp.depths = double(squeeze(temp.var_struct.grid.z(:, nearxy(temp.var_struct.grid.lon, temp.var_struct.grid.lat, temp.tspoints(1), temp.tspoints(2)))));
                            else
                                
                                temp.depths = temp.var_struct.grid.z;
                            end
                            
                      
                            temp.profile = double(squeeze(temp.var_struct.data(1, :, nearxy(temp.var_struct.grid.lon, temp.var_struct.grid.lat, temp.tspoints(1), temp.tspoints(2)))));
                            
                            figure;
                            plot(temp.profile, temp.depths)
                            xlabel(temp.var.attribute('units'))
                            ylabel('Depth')
                            title([temp.var.name, ' at ', num2str(temp.var_struct.grid.lon(nearxy(temp.var_struct.grid.lon, temp.var_struct.grid.lat, temp.tspoints(1), temp.tspoints(2)))), ', ', num2str(temp.var_struct.grid.lat(nearxy(temp.var_struct.grid.lon, temp.var_struct.grid.lat, temp.tspoints(1), temp.tspoints(2))))])
                            grid('on')
                        case 3
                            if isfield(temp.var_struct.grid, 'time')
                                [temp.temp temp.temp.temp temp.profile] = nearxy(temp.var_struct.grid.lon, temp.var_struct.grid.lat, temp.tspoints(1), temp.tspoints(2));
                                try
                                temp.data = double(squeeze(temp.var.data((temp.ind-100):(temp.ind+100), temp.profile(1), temp.profile(2))));
                                temp.depths = temp.var.grid_interop((temp.ind-100):(temp.ind+100), temp.profile(1), temp.profile(2));
                                catch
                                    try
                                        temp.data = double(squeeze(temp.var.data((temp.ind-100):end, temp.profile(1), temp.profile(2))));
                                        temp.depths = temp.var.grid_interop((temp.ind-100):end, temp.profile(1), temp.profile(2));
                                    catch
                                        try
                                            temp.data = double(squeeze(temp.var.data(1:(temp.ind+100), temp.profile(1), temp.profile(2))));
                                            temp.depths = temp.var.grid_interop(1:(temp.ind+100), temp.profile(1), temp.profile(2));
                                        catch
                                             temp.data = double(squeeze(temp.var.data(1:end, temp.profile(1), temp.profile(2))));
                                             temp.depths = temp.var.grid_interop(1:end, temp.profile(1), temp.profile(2));
                                        end
                                    end
                                end
                                figure;
                                plot(temp.depths.time, temp.data);
                                datetick('x','keeplimits');
%                                 xlabel(temp.var.attribute('units'))
                                ylabel([temp.var.name, ' in ', temp.var.attribute('units')]);
                                title([temp.var.name, ' at ', num2str(temp.var_struct.grid.lon(nearxy(temp.var_struct.grid.lon, temp.var_struct.grid.lat, temp.tspoints(1), temp.tspoints(2)))), ', ', num2str(temp.var_struct.grid.lat(nearxy(temp.var_struct.grid.lon, temp.var_struct.grid.lat, temp.tspoints(1), temp.tspoints(2))))])
%                                 grid('on')
                            else
                                % do nothing if z slice of 3d no time data chunk
                            end
                        case 2
                            % do nothing if 2d no time no z data surface
                    end
            case 'extend'
                out = get(gcf, 'UserData');
                out.data =[];
                out.depths = [];
                out.x = [];
                if ~isfield(out, 'path') || isempty(out.path)
                    warning('Please choose path by CTRL-LeftClicking on the model.');
                
                else
                [out.x, out.depths, out.data] = vsliceg(double(squeeze(out.var_struct.data)), out.var_struct.grid, out.path(:,1), out.path(:,2));
                out.path = [];
                set(gcf, 'UserData', out);
                figure;
                h = pcolor(out.x, out.depths, out.data);
                colorbar
                shading interp
                set(h, 'EdgeColor', 'k', 'EdgeAlpha', .1);
                title(out.var.attribute('units'))
                ylabel('Depth')
                xlabel('Distance in km');
                end

        end
        
    end
end

