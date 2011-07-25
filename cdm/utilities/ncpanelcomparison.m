% ncpanelcomparison - Function to compare 2 planer 2d horzitonal netcdf model/sat datasets.
function [] = ncpanelcomparison(dap1, dap2, colors, style, logscale)
if isempty(colors)
    colors = 'default';
end
if isempty(style)
    style = 'default';
end

        
switch nargin
    case 0
        datasets = ramadda_browse; 
    case 1
        error('Less than two arguments. Please enter two dataset opendap links. Or try using modellook.m.')
    otherwise
        switch isempty(dap1) || isempty(dap2)
            case 1
                datasets = ramadda_browse;
            case 0
                datasets.dataset_1 = ncgeodataset(dap1);
                datasets.dataset_2 = ncgeodataset(dap2);
        end
end

date = inputdlg('Input Time As YYYY-MM-DD-hh:mm:ss');

figure;
b = subplot(1,2,1);
panellook(datasets.dataset_1, date, colors, style);

% ar = get(b, 'PlotBoxAspectRatio');
% axis square
v = subplot(1,2,2);
panellook(datasets.dataset_2, date, colors, style);

% colorbar('PlotBoxAspectRatio', [1 40 1])
% set(v, 'PlotBoxAspectRatio', ar);
% axis square
linkaxes([b v]);
linkprop([b v], 'CLim');

%%
function [] = panellook(nc, date, colors, style)


% nc = ncgeodataset(dap);
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
out.var = var;
out.var_struct = var_struct;
datasets.colors = colors;
datasets.style = style;

set(gcf, 'UserData', datasets)


switch style
    case 'default'
        if length(var.size) > 3
            if abs(mean(var_struct.grid.z(1,~isnan(var_struct.grid.z(1,:,:))))) < abs(mean(var_struct.grid.z(end,~isnan(var_struct.grid.z(end,:,:)))));
                h = pcolor(squeeze(double(var_struct.grid.lon)), squeeze(double(var_struct.grid.lat)),...
                    squeeze(double(var_struct.data(1, 1, :,:))));
            else
                h = pcolor(squeeze(double(var_struct.grid.lon)), squeeze(double(var_struct.grid.lat)),...
                    squeeze(double(var_struct.data(1, end, :,:))));
            end
            title({nc.attribute('title');[var.name, ' in ', var.attribute('units'), ' on ', datestr(tind.time)]}, 'interpreter', 'none');
        elseif length(var.size) > 2
            if isfield(var_struct.grid, 'z')
                h = pcolor(squeeze(double(var_struct.grid.lon)), squeeze(double(var_struct.grid.lat)),...
                    squeeze(double(var_struct.data(1, :,:))));
                title({nc.attribute('title');...
                    [var.name, ' in ', var.attribute('units')]}, 'interpreter', 'none');
            elseif isfield(var_struct.grid, 'time')
                h = pcolor(squeeze(double(var_struct.grid.lon)), squeeze(double(var_struct.grid.lat)),...
                    squeeze(double(var_struct.data(1, :,:))));
                title({nc.attribute('title');[var.name, ' in ', var.attribute('units'),  'on ', datestr(tind.time)]}, 'interpreter', 'none');
            else
            end
        elseif length(var.size) > 1
            h = pcolor(squeeze(double(var_struct.grid.lon)), squeeze(double(var_struct.grid.lat)),...
                squeeze(double(var_struct.data(:,:))));
            title({nc.attribute('title');[var.name, ' in ', var.attribute('units')]}, 'interpreter', 'none');
        end
        shading interp
        set(h, 'EdgeColor', 'k', 'EdgeAlpha', .2);
    case 'contour'
        if length(var.size) > 3
            if abs(mean(var_struct.grid.z(1,~isnan(var_struct.grid.z(1,:,:))))) < abs(mean(var_struct.grid.z(end,~isnan(var_struct.grid.z(end,:,:)))));
                [C,h] = contourf(squeeze(double(var_struct.grid.lon)), squeeze(double(var_struct.grid.lat)),...
                    squeeze(double(var_struct.data(1, 1, :,:))));
            else
                [C,h] = contourf(squeeze(double(var_struct.grid.lon)), squeeze(double(var_struct.grid.lat)),...
                    squeeze(double(var_struct.data(1, end, :,:))));
            end
            title({nc.attribute('title');[var.name, ' in ', var.attribute('units'), ' on ', datestr(tind.time)]}, 'interpreter', 'none');
        elseif length(var.size) > 2
            if isfield(var_struct.grid, 'z')
                [C,h] = contourf(squeeze(double(var_struct.grid.lon)), squeeze(double(var_struct.grid.lat)),...
                    squeeze(double(var_struct.data(1, :,:))));
                title({nc.attribute('title');...
                    [var.name, ' in ', var.attribute('units')]}, 'interpreter', 'none');
            elseif isfield(var_struct.grid, 'time')
               [C,h] = contourf(squeeze(double(var_struct.grid.lon)), squeeze(double(var_struct.grid.lat)),...
                    squeeze(double(var_struct.data(1, :,:))));
                title({nc.attribute('title');[var.name, ' in ', var.attribute('units'),  'on ', datestr(tind.time)]}, 'interpreter', 'none');
            else
            end
        elseif length(var.size) > 1
            [C,h] = contourf(squeeze(double(var_struct.grid.lon)), squeeze(double(var_struct.grid.lat)),...
                squeeze(double(var_struct.data(:,:))));
            title({nc.attribute('title');[var.name, ' in ', var.attribute('units')]}, 'interpreter', 'none');
        end
end



set(h, 'UserData', out);

set(h, 'ButtonDownFcn', a)



colorbar('SouthOutside', 'PlotBoxAspectRatio', [40 1 1])
colormap(colors)


    function [out temp tspoint datasets] = press(src, evnt)
        sel_typ = get(gcbf,'SelectionType');
        switch sel_typ
            case 'alt'
                out = get(gcbo, 'UserData');
                temp = get(gca,'CurrentPoint');
                if ~isfield(out, 'path')
                    out.path = [];
                end
                out.path(end+1,1:2) = temp(1,1:2);
                set(gcbo, 'UserData', out);
            case 'open'
                tspoint = [];
                out = get(gca,'CurrentPoint');
                out = out(1,1:2);
                b = figure;
                for ll = 1:2;
                    figure(gcbf);
                    j = subplot(1,2,ll);
                    temp = get(j, 'Children');
                    temp = get(temp, 'UserData');
%                     temp.points = get(gca,'CurrentPoint');
%                     temp.tspoints = temp.points(1,1:2);
                    
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
                                
                                temp.depths = double(squeeze(temp.var_struct.grid.z(:, nearxy(temp.var_struct.grid.lon, temp.var_struct.grid.lat, out(1), out(2)))));
                            else
                                
                                temp.depths = temp.var_struct.grid.z;
                            end
                           
                            temp.profile = double(squeeze(temp.var_struct.data(1, :, nearxy(temp.var_struct.grid.lon, temp.var_struct.grid.lat, out(1), out(2)))));
                            
                            
                            figure(b);
                            v(ll) = subplot(1,2,ll);
                            if temp.depths(2)>0
                                plot(temp.profile, temp.depths.*-1)
                            else
                                plot(temp.profile, temp.depths)
                            end
                            
                            xlabel(temp.var.attribute('units'))
                            ylabel('Depth')
                            title([temp.var.name, ' at ', num2str(temp.var_struct.grid.lon(nearxy(temp.var_struct.grid.lon, temp.var_struct.grid.lat, out(1), out(2)))), ', ', num2str(temp.var_struct.grid.lat(nearxy(temp.var_struct.grid.lon, temp.var_struct.grid.lat, out(1), out(2))))])
                            grid('on')
                            tspoint = horzcat(tspoint, temp.profile);
                            
                        case 3
                            if isfield(temp.var_struct.grid, 'time')
                                temp.profile = [];
                                temp.data = [];
                                temp.depths = [];
                                [temp.temp temp.temp.temp temp.profile] = nearxy(temp.var_struct.grid.lon, temp.var_struct.grid.lat, out(1), out(2));
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
                                figure(b);
                                v(ll) = subplot(1,2,ll);
                                plot(temp.depths.time, temp.data);
                                datetick('x','keeplimits');
                                %                                 xlabel(temp.var.attribute('units'))
                                ylabel([temp.var.name, ' in ', temp.var.attribute('units')]);
                                title([temp.var.name, ' at ', num2str(temp.var_struct.grid.lon(nearxy(temp.var_struct.grid.lon, temp.var_struct.grid.lat, out(1), out(2)))), ', ', num2str(temp.var_struct.grid.lat(nearxy(temp.var_struct.grid.lon, temp.var_struct.grid.lat, out(1), out(2))))])
                                %                                 grid('on')
                            else
                                % do nothing if z slice of 3d no time data chunk
                                % add profile later, can stilldo profile as in case 4 above
                            end
                        case 2
                            % do nothing if 2d no time no z data surface
                    end
                    
                    
                end
                tspoint(isnan(tspoint)) = [];
                subplot(v(1));
                xlim([min(min(tspoint))-(std(tspoint)/2) max(max(tspoint))+(std(tspoint)/2)]);
                subplot(v(2));
                xlim([min(min(tspoint))-(std(tspoint)/2) max(max(tspoint))+(std(tspoint)/2)]);
                linkaxes(v);
                tspoint = [];
                tspoint = get(v, 'YLim');
                tspoint = horzcat(tspoint{1},[min(temp.depths(~isnan(temp.profile))) max(temp.depths(~isnan(temp.profile)))]);
                ylim([min(min(tspoint)) max(max(tspoint))]);
                    
                    clear temp
                    
                    
                    
            case 'extend'
                datasets = get(gcbf, 'UserData');
                out = get(gcbo, 'UserData');
              
                if ~isfield(out, 'path') || isempty(out.path)
                    warning('Please choose path by CTRL-LeftClicking on the model.');
                
                else
                
                
                temp = out.path; 
                
                out.path = [];
                set(gcbo, 'UserData', out);
                
                figure(gcbf);
                j = subplot(121);
                out = get(j, 'Children');
                out = get(out, 'UserData');
                out.data =[];
                out.depths = [];
                out.x = [];
                [out.x, out.depths, out.data] = vsliceg(double(squeeze(out.var_struct.data)), out.var_struct.grid, temp(:,1), temp(:,2));
                
                q = figure;
                b = subplot(2,1,1);
                switch datasets.style
                    case 'contour'
                        if mean(out.depths) > 0;
                            [C,h] = contourf(out.x, out.depths.*-1, out.data);
                            tspoint = out.depths(~isnan(out.data)).*-1;
                        else
                            [C, h] = contourf(out.x, out.depths, out.data);
                            tspoint = out.depths(~isnan(out.data));
                        end
                    otherwise
                        if out.depths(2,2) > 0;
                            h = pcolor(out.x, out.depths.*-1, out.data);
                            tspoint = out.depths(~isnan(out.data)).*-1;
                        else
                            h = pcolor(out.x, out.depths, out.data);
                            tspoint = out.depths(~isnan(out.data));
                        end
                        shading interp
                        set(h, 'EdgeColor', 'k', 'EdgeAlpha', .1);
                end
                colorbar('PlotBoxAspectRatio', [1 45 1]);
                colormap(datasets.colors)
    
                
                title([out.var.attribute('units'), ' - ', out.var.dataset.attribute('title')])
                ylabel('Depth')
                
                
                figure(gcbf);
                j = subplot(122);
                out = get(j, 'Children');
                out = get(out, 'UserData');
                
                out.data =[];
                out.depths = [];
                out.x = [];
                
                [out.x, out.depths, out.data] = vsliceg(double(squeeze(out.var_struct.data)), out.var_struct.grid, temp(:,1), temp(:,2));
                
                figure(q)
                v = subplot(2,1,2);
                
                switch datasets.style
                    case 'contour'
                        if mean(out.depths) > 0;
                            [C,h] = contourf(out.x, out.depths.*-1, out.data);
                            tspoint = vertcat(reshape(tspoint, [numel(tspoint), 1]), reshape(out.depths(~isnan(out.data)).*-1,[numel(out.depths(~isnan(out.data))), 1]));
                        else
                            [C,h] = contourf(out.x, out.depths, out.data);
                            tspoint = vertcat(reshape(tspoint, [numel(tspoint), 1]), reshape(out.depths(~isnan(out.data)),[numel(out.depths(~isnan(out.data))), 1]));
                        end
                    otherwise
                        if out.depths(2,2) > 0;
                            h = pcolor(out.x, out.depths.*-1, out.data);
                            tspoint = vertcat(reshape(tspoint, [numel(tspoint), 1]), reshape(out.depths(~isnan(out.data)).*-1,[numel(out.depths(~isnan(out.data))), 1]));
                        else
                            h = pcolor(out.x, out.depths, out.data);
                            tspoint = vertcat(reshape(tspoint, [numel(tspoint), 1]), reshape(out.depths(~isnan(out.data)),[numel(out.depths(~isnan(out.data))), 1]));
                        end
                        shading interp
                        set(h, 'EdgeColor', 'k', 'EdgeAlpha', .1);
                end
                
                colorbar('PlotBoxAspectRatio', [1 45 1]);
                colormap(datasets.colors)
                
                
                
                title([out.var.attribute('units'), ' - ', out.var.dataset.attribute('title')])
                ylabel('Depth')
                
                subplot(b)
                ylim([min(min(tspoint))-1 max(max(tspoint))]);
                temp = [];
                temp(1,:) = get(b, 'CLim');
                subplot(v)
                ylim([min(min(tspoint))-1 max(max(tspoint))]);
                xlabel('Distance in km');
                temp(2,:) = get(v, 'CLim');
                linkaxes([b v]);
                linkprop([b v], 'CLim');
                caxis([min(min(temp)) max(max(temp))])
                end

        end
        
    end
end

end