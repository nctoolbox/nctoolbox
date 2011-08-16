function [tspoint] = sectionpanel(style, out)
%
%
%
%
%
switch style
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

end