% MODELLOOK - Infrastructure for quick model visualization.
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

[select, ok] = listdlg('ListString', sizes,'SelectionMode','sing le', 'Name', 'ModelLook: Quick nc vis',...
    'PromptString', 'Please Select Variable:' );
var = nc.geovariable(vars{select});
tind = var.timewindow(date);
var_struct.data = var.data(tind.index, :, :, :);
var_struct.grid = var.grid_interop(tind.index, :, :, :);

a = @press;
figure;
out.var = var;
out.var_struct = var_struct;
out.ind = tind.index;

set(gcf, 'UserData', out);
if var.size > 3
    h = pcolor(squeeze(double(var_struct.grid.lon)), squeeze(double(var_struct.grid.lat)),...
        squeeze(double(var_struct.data(1, 1, :,:))));
elseif var.size > 2
    h = pcolor(squeeze(double(var_struct.grid.lon)), squeeze(double(var_struct.grid.lat)),...
        squeeze(double(var_struct.data(1, :,:))));
end
set(h, 'ButtonDownFcn', a)
% shading flat
colorbar
title({'ModelLook:';'Ctl-Click to create a path for vertical section and Shift-Click to plot.';...
    'Double-Click to create a profile plot at model location.';...
    ' '; [var.name, ' in ', var.attribute('units'), ' on ', datestr(tind.time)]});

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
%                 for out = 1:10:temp.size(1)
                    temp.data = [];
                    temp.depths = [];
                    temp.profile = [];
%                     try
%                         temp.data = temp.var.data(temp.ind, :, :, :);
                        temp.profile = interptoxy(double(squeeze(temp.var_struct.data)), temp.var_struct.grid.lon,...
                            temp.var_struct.grid.lat, temp.tspoints(1), temp.tspoints(2), 'linear');
                        if length(size(temp.var_struct.grid.z)) > 2
                            temp.depths = interptoxy(double(squeeze(temp.var_struct.grid.z)),...
                                temp.var_struct.grid.lon, temp.var_struct.grid.lat, temp.tspoints(1),...
                                temp.tspoints(2), 'linear');
                        else 
                            temp.depths = temp.var_struct.grid.z;
                        end
%                     catch
%                         temp.data = temp.var.data(out:end, :, :, :);
%                         temp.ts(:,out:end) = interptoxy(temp.data, temp.var_struct.grid.lon,...
%                             temp.var_struct.grid.lat, temp.points(1), temp.points(2), 'nearest');
%                     end
%                     temp.time = temp.var.
%                 end
                figure;
                plot(temp.profile, temp.depths)
                xlabel(temp.var.attribute('units'))
                ylabel('Depth')
                title([temp.var.name, ' at ', num2str(temp.tspoints(1)), ', ', num2str(temp.tspoints(2))])
                grid('on')
%                 for out = 1:length(temp.ts(:,1))
%                     plot(time.ts(out, :));
%                     hold on
%                 end
%                 hold off
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
                pcolor(out.x, out.depths, out.data)
                colorbar
                title(out.var.attribute('units'))
                ylabel('Depth')
                end

        end
        
    end
end

