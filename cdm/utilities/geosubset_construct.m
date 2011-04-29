function [ s ] = geosubset_construct(starttime, stoptime, t_stride, zmin, zmax, z_stride, eastmin, ...
  eastmax, x_stride, northmin, northmax, y_stride)
% SUBSET_STRUCT this is a function that, in a matlab way, creates a subset structure for nctoolbox's 
% geosubset or timegeosubset functions.
% >> geosubset_construct(1, 1, [], 1, 1, [], -77, -75.5, [],  36, 38, 1)
%
% Note: the starttime and stoptime arugments also accept matlab time formats datevec and datenum.
%

% Bounding Dimensions
s.lat = [northmin northmax]; % both ends must be explicit
s.lon = [eastmin eastmax]; % both ends must be explicit
s.z_index = [zmin zmax]; % both ends must be explicit
s.time = {starttime stoptime}; %both ends must be explicit

% Stride Arugments
if isempty(z_stride) % Z
  s.z_stride = 1;
else
  s.z_stride = z_stride;
end

if isempty(t_stride) % Time
  s.t_stride = 1;
else
  s.t_stride = t_stride;
end

if isempty(x_stride) %lon 
  s.xy_stride(1) = 1;
else
  s.xy_stride(1) = x_stride;
end

if isempty(y_stride) %lat
  s.xy_stride(2) = 1;
else
  s.xy_stride(2) = y_stride;
end


end

