function [hax,hax_noxlabel,hax_noylabel,hax_left,hax_bottom] = nfigaxes(naxes,gaps,xlim,ylim)
% Create a set of axes on a single page
% [hax,hax_noxlabel,hax_noylabel,hax_left,hax_bottom] = ...
%                                         nfigaxes(naxes,gaps,xlim,ylim)
%     Inputs:
%           NAXES is a vector setting number of axes to create [ny nx] 
%     Optional inputs:  If missing or empty [] the defaults are used
%           GAPS [xgap ygap] overrides the default gaps between axes
%           XLIM [left right] overrides the default limits of the page 
%           positions within which the axes lie
%           YLIM [bot top] overrides the default limits of the page 
%           positions within which the axes lie
%           (XLIM,YLIM) must be within 0->1
%           defaults:
%           left and right x limits, top and bottom y limits
%           xl = 0.12; xr = 0.88; yb = 0.2; yt = 0.8;
%           gaps between axes
%           xgap = 0.02; ygap = 0.05;
%
%     Output:
%           HAX are the handles for the axes
%           HAX_NOXLABEL,HAX_NOYLABEL are the handles of the interior 
%           axes that may clutter the plot, and the labels can be removed
%           with
%           set(hax_noxlabel,'xtick',[])
%           set(hax_noylabel,'ytick',[])
%           HAX_LEFT,HAX_RIGHT are the handles of the axes you may wish to
%           modify the labels on
%
%     The axes number from top left corner going across x dir first
%     e.g.  1  -> 2 -> 3
%           4  -> 5 -> 6
%     Use command axes(hax(i)) to select axes in which to plot
%
% John Wilkin Oct 96

% defaults:
% left and right x limits, top and bottom y limits
xl = 0.12; xr = 0.88; yb = 0.2; yt = 0.8;
% gaps between axes
xgap = 0.02; ygap = 0.05;

if nargin > 1
  if isstr(gaps)
    if strcmp(gaps,'max')
      xl = 0.05; xr = 0.95; yb = 0.05; yt = 0.95;
      gaps = [0 0];
    elseif strcmp(gaps,'full')
      xl = 0; xr = 1; yb = 0; yt = 1;
      gaps = [0 0];      
    end
  end
  % gap between axes
  if ~isempty(gaps)
    xgap = gaps(1);
    ygap = gaps(2);
  end
  if nargin > 2
    if ~isempty(xlim)
      xl = xlim(1);
      xr = xlim(2);
    end
    if nargin > 3
      if ~isempty(ylim)
	yb = ylim(1);
	yt = ylim(2);
      end
    end
  end
end

% compute size of axes to fit these limits
dx = ((xr-xl)-(naxes(2)-1)*xgap)/naxes(2);
dy = ((yt-yb)-(naxes(1)-1)*ygap)/naxes(1);

% build matrix of position vectors
pos = NaN*ones([prod(naxes) 4]);
nax = 0;
for yax=naxes(1)-1:-1:0
  for xax=0:naxes(2)-1
    nax = nax+1;
    pos(nax,1) = xl + xax*(dx+xgap);
    pos(nax,2) = yb + yax*(dy+ygap);
  end
end
pos(:,3) = dx*ones([prod(naxes) 1]);
pos(:,4) = dy*ones([prod(naxes) 1]);

% create axes and save handles
hax = [];
for i=1:prod(naxes)
  hax(i)=axes('position',pos(i,:));
end

if nargout > 1
  % separate the handles of the axes that should have the labels
  % removed to avoid clutter
  hax_bottom = (prod(naxes)-naxes(2)+1):prod(naxes);
  hax_noxlabel = hax;
  hax_noxlabel(hax_bottom) = [];
  
  hax_left = 1:naxes(2):prod(naxes);
  hax_noylabel = hax;
  hax_noylabel(hax_left) = [];

  hax_left = hax(hax_left);
  hax_bottom = hax(hax_bottom);
  
end
