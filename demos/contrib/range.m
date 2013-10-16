function y = range(x)
%RANGE  Calculates the difference between the maximum and minimum values. 
% 
% CALL: r = range(x) 
%
%   r =  max(x) - min(x), i.e. a scalar or vector giving the range of x
%   x = vector or matrix
%
% Example:
%   x=wexprnd(5,1,10)
%   [min(x) max(x) range(x)]
%
% See also: max, min, iqr

% By pab 06.02.2000

%y = max(x) - min(x);
y = [min(x(:)) max(x(:))]
