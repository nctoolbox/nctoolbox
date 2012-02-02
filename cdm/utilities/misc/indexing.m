function [start last stride] = indexing(indices, varShape)

%mVariable/parseIndices - parses subscripting references into
% start,last and stride based on variable shape.

% [start last stride] = parseIndices(indices, varShape)
% where,
% indices - array subscripted input to index type '()'. e,g (:,1:2,end)
% varShape - variable shape

% Sachin Kumar Bhate (skbhate@ngi.msstate.edu) (C) 2008
% Mississippi State University
%
% Note: Logic to calculate start/last/stride partially borrowed
% from ncvar/subsref ['netcdf-toolbox' by Chuck Dhenam (c) 1996-97].
% and modified to suit our needs.
%
% Modified by Alexander Crosby 2010, 2011 (nctoolbox)
% http://code.google.com/p/nctoolbox

%init
start = ones(1, length(varShape));
last = ones(1, length(varShape)); 
stride = ones(1, length(varShape));

if nargin < 2, help(mfilename), return, end

try
    for i = 1:length(indices)
        if isa(indices{i}, 'double')
            if any(diff(diff(indices{i})))
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
                %last(i) = length(k);
                d = 0;
                if length(k) > 1, d = diff(k); end
                stride(i) = max(d(1), 1);
                last(i) = k(end);%*stride(i);
                %maxlast1 = fix((varShape(i)-start(i)+stride(i)) ./ stride(i));
            else
                last(i) = -1;
                if i == length(indices) && i < length(varShape)
                    j = i+1:length(varShape);
                    last(j) = -ones(1, length(j));
                end
            end
        end
        %         start(start < 1) = 1;
        stride(stride < 1) = 1;
        for i = 1:length(last)
            if last(i) < 0
                %             if last(i) == -1
                %                 maxlast = fix((varShape(i)-start(i)+stride(i)) ./ stride(i)); % take a look at this
                maxlength = varShape(i);
                last(i) = maxlength + last(i) + 1;
            end
            if start(i) < 0
                %             if start(i) == -1
                %                 maxlast = fix((varShape(i)-start(i)+stride(i)) ./ stride(i)); % take a look at this
                maxlength = varShape(i);
                start(i) = maxlength + start(i) + 1;
            end
            
        end
    end
catch %gets the last error generated
    err = lasterror();
    disp(err.message);
end
end
