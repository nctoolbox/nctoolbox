function v = value4key(a, key)
% VALUE4KEY - Lookup a value from a 2-by-n cell array 
%
% Use as: 
%   v = value4key(a, key)
%
% Inputs:
%   a   = a 2-by-n cell-array of strings. (e.g. as returned by
%         ncdataset.attributes). The first column should represent the
%         'keys', the second column contains the corresponding 'value'.
%   key = The key in the first column to search for. For example,
%         if you want the value for the 'units' attribute, the key would be
%         'units'
%
% Output:
%   v = The value for the corresponding key. If no match is found an
%       empty char array is returned.
%
% Example:
%   ds = ncdataset('http://dods.mbari.org/cgi-bin/nph-nc/data/ssdsdata/deployments/m1/200810/OS_M1_20081008_TS.nc');
%   at = ds.attributes('TEMP');
%   units = value4key(at, 'units'); % returns 'celsius'

keys = a(:, 1);    % Cell array of keys
values = a(:, 2);  % Cell array of values
mtch = ismember(keys, key);
if any(mtch)
    v = values{mtch};  % Retrieve the units value
else
    v = '';
end
