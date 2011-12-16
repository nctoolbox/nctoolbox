function [h]=arrows(x,y,w,fac,color)
% function [h]=arrows(x,y,w,fac,color)
%    Draws arrows with their tails at each point corresponding
%	to identical indices in the matrices x,y.  The matrices u
%	and v are the components of the vector to be represented.
%	fac is the scaling factor 
%
%  Geometry of arrowheads (choosing HEADA and HEADL):
%    If the arrow is defined by the points A B C B D where A is the base of 
%  the arrow, B is the head, and C and D are the corners of the arrowhead, then
%  HEADA is the angle BAC (or BAD), and HEADL is the ratio of distances AC/AB.
%

% revision 3/20/97 to use nans for line breaks
% much more efficient and returns only a single handle
  
HEADA=10*pi/180; HEADL=.75;
z=x(:)+i*y(:);

if nargin < 5,color='red';end
if nargin < 4,help arrows,end
w=w(:)*fac;
r=w*HEADL; wr1=r*exp(+i*HEADA); wr2=r*exp(-i*HEADA);
wplot=ones(length(z),6); 
wplot(:,1)=z; 
wplot(:,[2,4])=(z+w)*ones(1,2);
wplot(:,3)=z+wr1; wplot(:,5)=z+wr2;
wplot(:,6)=z*nan;
wplot=wplot.';
wplot=wplot(:);
%z=eps*ones(size(wplot));
%h=line(real(wplot),imag(wplot),z,'color',color);
h=line(real(wplot),imag(wplot),'color',color);
set(h(1),'userdata',fac);
