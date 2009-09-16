function sdn = utc2sdn(utc)
% UTC2SDN   - Convert UTC seconds to matlab date format
%
% Use as: sdn = utc2sdn(utc);
%
% sdn = Matlab datenumber (se DATENUM and DATESTR)
% utc = Seconds since 01 Jan 1970 00:00:00

% Brian Schlining
% 12 Apr 2000

%datenum('01 Jan 1970 00:00:00') = 719529
sdn = utc/60/60/24 + 719529;