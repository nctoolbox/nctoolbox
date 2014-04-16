function [start, last, stride] = indexing(indices, varShape)

% INDEXING - parses subscripting references into
% start,last and stride based on variable shape.
%
% Use as:
%   [start last stride] = indexing(indices, varShape)
% 
% Arguments:
%   indices - Array subscripted input to index type '()'. e,g (:,1:2,end)
%   varShape - Variable shape
%
% Return:
%   start - Starting index. Often called 'first'
%   last - Ending index.
%   stride = stride for each dimension of the array
% 
%   start, last and stride will all be the same size as 'varShape'


% Sachin Kumar Bhate (skbhate@ngi.msstate.edu) (C) 2008
% Mississippi State University
%
% Note: Logic to calculate start/last/stride partially borrowed
% from ncvar/subsref ['netcdf-toolbox' by Chuck Denham (c) 1996-97].
% and modified to suit our needs.
%
% Modified by Alexander Crosby 2010, 2011 (nctoolbox)
% https://github.com/nctoolbox/nctoolbox

%init
start = ones(1, length(varShape));
last = ones(1, length(varShape));
stride = ones(1, length(varShape));

if nargin < 2, help(mfilename), return, end

try
    for i = 1:length(indices)
        if isa(indices{i}, 'double')
            if any(diff(diff(indices{i})))
                %indices{i}
                error('MATLAB:mVariable:parseIndices', ...
                    'Indexing strides must be positive and constant.')
            end
        end
    end
    
    if prod(varShape) > 0
        for i = 1:min(length(indices), length(varShape))
            k = indices{i};
            if ~ischar(k) && ~strcmp(k, ':') && ~strcmp(k, '-')
                start(i) = k(1);
                d = 0;
                if length(k) > 1
                    d = diff(k); 
                end
                stride(i) = max(d(1), 1);
                last(i) = k(end);
            else
                last(i) = -1;
                if i == length(indices) && i < length(varShape)
                    j = i + 1:length(varShape);
                    last(j) = -ones(1, length(j));
                end
            end
        end
        stride(stride < 1) = 1;
        for i = 1:length(last)
            if last(i) < 0
                maxlength = varShape(i);
                last(i) = maxlength + last(i) + 1;
            end
            if start(i) < 0
                maxlength = varShape(i);
                start(i) = maxlength + start(i) + 1;
            end
            
        end
    end
catch me
    % gets the last error generated
    fprintf(1, '%s occured on the following indices:', me.message);
    disp(indices)
    error('NCTOOLBOX:indexing', me.message);
end

