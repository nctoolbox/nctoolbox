function i=date_index(jdmat,dn1,dn2);
% date_index - Finds indices of DATENUM array in specified date range
% Usage: i=date_index(jdmat,dn1,[dn2]);
%   Inputs: jdmat = Matlab datenum array
%            dn1 = Matlab datenum, datestr or datevec
%            dn2 = Matlab datenum, datestr or datevec (optional)
%   Output: index of closest date to dn1 (if dn2 is not specified)
%           indices between dn1 and dn2 (if dn2 is specified) (returns NaN
%           if no dates are in range)
%   Examples:
%    i=date_index(jd,[1990 4 5 0 0 0]);
%    i=date_index(jd,'1990-Apr-5 00:00');
%    returns the index of jd that is closest to to April 5, 1990 00:00
%
%    ii=date_index(jd,[1990 4 5 0 0 0],[1992 3 1 12 0 0]);
%    returns the indices of jd that are between 0000 April 5, 1990
%      and 1200 March 1, 1992.
%
%  rsignell@usgs.gov

i=NaN;
jdmat=jdmat(:);
dn1=datenum(dn1);
% if(dn1<jdmat(1) | dn1>jdmat(end))
%     sprintf('Error: Requested date %s is outside the range %s - %s.',...
%         datestr(jdmat(1)),datestr(jdmat(end)),datestr(dn1));
%     return
% end
switch nargin
    case 2
        i=near(jdmat,dn1,1);
    case 3
        dn2=datenum(dn2);
%         if(dn2<jdmat(1) | dn2>jdmat(end)),
%             sprintf('Error: Requested date %s is outside the range %s - %s.',...
%                 datestr(jdmat(1)),datestr(jdmat(end)),datestr(dn2));
%             return
%         end
        i=find(jdmat>=dn1&jdmat<=dn2);
    otherwise
        help date_index
end

