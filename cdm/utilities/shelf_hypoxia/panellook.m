%%
function [modeldata] = panellook(nc, var, date, colors, style, units, varargin)

if ~isempty(varargin)
    scale = varargin{2};
    depth = varargin{1};
else
end

var = nc.geovariable(var);

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

% Units conversion dance:
% Try to convert to the required input units:
try
    var_struct.data = ncunits(var_struct.data, var.attribute('units'), units);
catch
    var_struct.data = hypoxia_units(var_struct.data, var.attribute('units'), units);
end


a = @press;
out.var = var;
out.var_struct = var_struct;
out.units = units;
datasets.colors = colors;
datasets.style = style;
modeldata = out;

set(gcf, 'UserData', datasets)


switch style
    case 'default'
        h = default_panel(nc, var, var_struct, tind, depth, units);
    case 'contour'
        h = contour_panel(nc, var, var_struct, tind, depth, units);
end

% Explicitly add model data to axes and set the click/button-down functionality to the axes
set(h, 'UserData', out);
set(h, 'ButtonDownFcn', a);

% Set location and shape of colorbar
colorbar('SouthOutside', 'PlotBoxAspectRatio', [40 1 1])
colormap(colors) % Set colormap to use
if ~isempty(varargin{1})
    switch scale
        case 'linear'
            % do nothing
        case 'log'
            colorbar_log([10e-1 10e2]);
    end
end


    function [out temp tspoint datasets] = press(src, evnt)
        % This function is used to handle the button down events on the axes
        % Because of the way the variables are created and handled in the function handle
        % and reused through the axes' user data many of the variables are used to represent fields they 
        % are not named to represent.
        %
        %
        %
        sel_typ = get(gcbf,'SelectionType');
        switch sel_typ
            case 'alt' % Right Click
                out = get(gcbo, 'UserData');
                temp = get(gca,'CurrentPoint');
                if ~isfield(out, 'path')
                    out.path = [];
                end
                out.path(end+1,1:2) = temp(1,1:2);
                set(gcbo, 'UserData', out);
            case 'open' % Double click
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
                             % Get profile at point from 4d (t,z,y,x) variable
                             temp = profilefrom4d(temp, out);                   
                            
                             % Plot
                            figure(b);
                            v(ll) = subplot(1,2,ll);
                            if temp.depths(2)>0 % Double check orientation of vertical coords
                                plot(temp.profile, temp.depths.*-1)
                            else
                                plot(temp.profile, temp.depths)
                            end
                            
                            % Add labels 
                            xlabel(temp.units)
                            ylabel('Depth')
                            title([temp.var.name, ' at ', num2str(temp.var_struct.grid.lon(nearxy(temp.var_struct.grid.lon, temp.var_struct.grid.lat, out(1), out(2)))), ', ', num2str(temp.var_struct.grid.lat(nearxy(temp.var_struct.grid.lon, temp.var_struct.grid.lat, out(1), out(2))))])
                            grid('on')
                            tspoint = horzcat(tspoint, temp.profile); % temp variable to calc plot extents
                             
                        case 3
                            if isfield(temp.var_struct.grid, 'time')
                                temp.profile = [];
                                temp.data = [];
                                temp.depths = [];
                                
                                %Get time series from a 3d variable (t, y, x)
                                temp = timeseriesfrom3d(temp, out);
                                
                                % Plot
                                figure(b);
                                v(ll) = subplot(1,2,ll);
                                plot(temp.depths.time, temp.data);
                                datetick('x','keeplimits');
                                %                                 xlabel(temp.var.attribute('units'))
                                ylabel([temp.var.name, ' in ', temp.units]);
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
                tspoint(isnan(tspoint)) = []; % Remove nan
                % Set x scale/extents
                subplot(v(1));
                xlim([min(min(tspoint))-(std(tspoint)/2) max(max(tspoint))+(std(tspoint)/2)]);
                subplot(v(2));
                xlim([min(min(tspoint))-(std(tspoint)/2) max(max(tspoint))+(std(tspoint)/2)]);
                
                linkaxes(v);% Link both axes for zooming and panning
                
                tspoint = []; % maybe i dont user the tspoint var above, i should check
                tspoint = get(v, 'YLim');
                tspoint = horzcat(tspoint{1},[min(temp.depths(~isnan(temp.profile))) max(temp.depths(~isnan(temp.profile)))]);
                ylim([min(min(tspoint)) max(max(tspoint))]); % set the y scale/extents of the profile plots
                    
                clear temp
                    
                    
                    
            case 'extend' % Center click
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
                % Create data for section
                [out.x, out.depths, out.data] = vsliceg(double(squeeze(out.var_struct.data)), out.var_struct.grid, temp(:,1), temp(:,2), 'linear');
                
                % Plot
                q = figure;
                b = subplot(2,1,1);
                tspoint = sectionpanel(datasets.style, out); %plot the section with a style
                colorbar('PlotBoxAspectRatio', [1 45 1]);
                colormap(datasets.colors)
                title([out.var.name, ' in ',out.units,  ' - ', out.var.dataset.attribute('title')])
                ylabel('Depth')
                
                
                figure(gcbf);
                j = subplot(122);
                out = get(j, 'Children');
                out = get(out, 'UserData');
                
                out.data =[];
                out.depths = [];
                out.x = [];
                
               % Create data for section
                [out.x, out.depths, out.data] = vsliceg(double(squeeze(out.var_struct.data)), out.var_struct.grid, temp(:,1), temp(:,2), 'linear');
                
                % Plot
                figure(q)
                v = subplot(2,1,2); 
                tspoint = sectionpanel(datasets.style, out); %plot the section with a style
                colorbar('PlotBoxAspectRatio', [1 45 1]);
                colormap(datasets.colors)
                title([out.var.name, ' in ',out.units, ' - ', out.var.dataset.attribute('title')])
                ylabel('Depth')
                
                % Set scale/extents, link axis, add labels
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