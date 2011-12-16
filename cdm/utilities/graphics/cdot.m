function cdot(x,y,z,c,dot_size,label,range);
% CDOT plots colored dots over the range of data
% 
% Usage:  cdot(x,y,z,c,dot_size,label,range);
%
%   x=vector of x locations
%   y = vector of x locations
%   z = vector of z values
%   c = colormap
%   dot_size = Marker size for dots (e.g 8)
%   label = 1 to label with value, 0 otherwise
%   range = range of data for colormap
%
%  example 
%    [x,y,z]=peaks(20);cdot(x(:),y(:),z(:),jet,20,1,[-8 8])
%
 [m,n]=size(c);
 ndots=length(z);
 zsave=z;
 if(exist('range')==1),
  zmin=range(1);
  zmax=range(length(range))+eps;
  iz=find(z>=zmax);
  if(~isempty(iz)),
    z(iz)=(zmax-eps)*ones(size(iz));
  end
  iz=find(z<zmin);
  if(~isempty(iz)),
    z(iz)=zmin*ones(size(iz));
  end
 else
  zmin=min(z(:));
  zmax=max(z(:))+eps;
 end
 zinc=(zmax-zmin)/m;  
% set(gca,'xlim',[min(x) max(x)],'ylim',[min(y) max(y)]);
 for i=1:m;
   z1=zmin+zinc*(i-1);
   z2=zmin+zinc*i;
   ind=find(z>=z1& z<=z2);
   line(x(ind),y(ind),'linestyle','none','marker','.','markersize',...
      dot_size,'color',c(i,:));
 end
 if(label==1);
     for i=1:ndots
%      text(x(i),y(i),sprintf(' %d',zsave(i)),'fontsize',14,...
      text(x(i),y(i),sprintf('%5.2f',zsave(i)),...
      'HorizontalAlignment','left','VerticalAlignment','bottom')
%      bgtext(x(i),y(i),sprintf(' %d',zsave(i)),[1 1 .9],'fontsize',10,...
%       bgtext(x(i),y(i),sprintf(' %d',zsave(i)),[1 1 .9]);
    end
 end
