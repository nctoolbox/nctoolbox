% ARG2HASH - Varargin to 2 by x cell hash
function [cellout] = arg2hash(in)


    c = 1;
    for i = 1:2:length(in)
        cellout{c,1} = in{i};
        cellout{c,2} = in{i+1};
        
        c = c+1;
    end








end