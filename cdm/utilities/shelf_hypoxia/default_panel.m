function [h] = default_panel(nc, var, var_struct, tind, depth, units)
% DEFAULT_PANEL
%
%
%

if length(var.size) > 3
    switch depth
        case 'surface'
            if abs(mean(var_struct.grid.z(1,~isnan(var_struct.grid.z(1,:,:))))) < abs(mean(var_struct.grid.z(end,~isnan(var_struct.grid.z(end,:,:)))));
                h = pcolor(squeeze(double(var_struct.grid.lon)), squeeze(double(var_struct.grid.lat)),...
                    squeeze(double(var_struct.data(1, 1, :,:))));
            else
                h = pcolor(squeeze(double(var_struct.grid.lon)), squeeze(double(var_struct.grid.lat)),...
                    squeeze(double(var_struct.data(1, end, :,:))));
            end
            title({nc.attribute('title');[var.name, ' in ', units, ' on ', datestr(tind.time)]}, 'interpreter', 'none');
        otherwise %use zsliceg to interpolate to z value
            slice = zsliceg(squeeze(double(var_struct.data(1, :, :,:))), var_struct.grid.z, depth);
            h = pcolor(squeeze(double(var_struct.grid.lon)), squeeze(double(var_struct.grid.lat)),...
                    slice);
            title({nc.attribute('title');[var.name, ' in ', units, ' on ', datestr(tind.time),  ' at ', num2str(depth)]}, 'interpreter', 'none');
    end
            
    elseif length(var.size) > 2
    if isfield(var_struct.grid, 'z')
        h = pcolor(squeeze(double(var_struct.grid.lon)), squeeze(double(var_struct.grid.lat)),...
            squeeze(double(var_struct.data(1, :,:))));
        title({nc.attribute('title');...
            [var.name, ' in ', units]}, 'interpreter', 'none');
    elseif isfield(var_struct.grid, 'time')
        h = pcolor(squeeze(double(var_struct.grid.lon)), squeeze(double(var_struct.grid.lat)),...
            squeeze(double(var_struct.data(1, :,:))));
        title({nc.attribute('title');[var.name, ' in ', units,  ' on ', datestr(tind.time)]}, 'interpreter', 'none');
    else
    end
elseif length(var.size) > 1
    h = pcolor(squeeze(double(var_struct.grid.lon)), squeeze(double(var_struct.grid.lat)),...
        squeeze(double(var_struct.data(:,:))));
    title({nc.attribute('title');[var.name, ' in ', units]}, 'interpreter', 'none');
end
        shading interp
        set(h, 'EdgeColor', 'k', 'EdgeAlpha', .2);
        
end