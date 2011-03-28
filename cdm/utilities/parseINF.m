function [start,count, stride]=parseINF(start, count, stride, shape)

%Utilties/parseINF - Assign 'inf' the correct dimension value 
% based on variable shape.

% [start, count, stride] = parseINF(start, count, stride, shape)
% where,
%           start = A 1-based array specifying the position in the file to begin reading or corner of a
%                   hyperslab e.g [1 2 1]. Specify 'inf' for last index.
%                   The values specified must not exceed the size of any
%                   dimension of the data set.
%           count = A 1-based array specifying the length of each dimension to read. e.g [6 10 inf]. 
%                   Specify 'inf' to get all the values from start index.
%           stride = A 1-based array specifying the interval between the
%                    values to read. e.g [1 2 2] 
% shape - variable shape

% Sachin Kumar Bhate (skbhate@ngi.msstate.edu) (C) 2008
% Mississippi State University
%


    if nargin < 4, help(mfilename), return, end

    try 
        startInf = find(isinf(start));
        countInf = find(isinf(count));

        if (~isempty(startInf)) %contains inf
            for i=1:length(startInf)
                start(startInf(i)) = shape(startInf(i));
                count(startInf(i)) = 1;
            end
        end
        if (~isempty(countInf)) %contains inf
            for i=1:length(countInf)
                count(countInf(i)) = shape(countInf(i))- start(countInf(i)) + 1;
            end 
        end
        
        %fix the final count based on the stride
        newcount = ones(1, length(shape));
        for i=1:length(shape)
            maxcount = fix((shape(i)-start(i)+stride(i)) ./ stride(i));          
            newcount(i)= fix((count(i)-1)/stride(i)) + 1;                
            if (newcount(i) > maxcount)                
                error('MATLAB:Utilites:parseINF',...
                        'Invalid count "%d", exceeds dimension "%d"', newcount(i),maxcount); 
            end
        end
        count = newcount;
        
    catch %gets the last error generated 
        err = lasterror();    
        disp(err.message);
    end 
end
