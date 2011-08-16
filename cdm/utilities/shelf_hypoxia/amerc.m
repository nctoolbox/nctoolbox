function [xrng,yrng] = amerc
% set DataAspectRatio to Mercator proportions at figure centre
ylim = get(gca,'ylim');
set(gca,'DataAspectRatio',[1 cos(mean(ylim)*pi/180) 1]);
if nargout~=0
    % output box dimensions in km
    yrng = sw_dist(ylim,[0 0],'km');
    xrng = sw_dist(mean(get(gca,'ylim'))*[1 1],get(gca,'xlim'),'km');
end
