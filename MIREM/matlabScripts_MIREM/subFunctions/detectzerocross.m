function [zerocross, set]  = detectzerocross(x)
    % computes the zerocrossing indexes and if the crossing is an
    % onset/offset. 
    % Values returned :
    % zerocross : a vector of the indexes in the vector x right before the
    % zero-crossing happens.
    % 
    % set : return a vector the same length as zerocross. This vector is
    % only composed of 1 and -1. It is related to the zerocross vector, 1
    % means this crossing is an onset crossing and -1 means it's an offset
    % crossing


    idx                    = [1:length(x)];
    last                   = circshift(x,1);
    zerocross              = idx(x.*last<0);  %  = TRUE when x crosses zero
    set                    = ones(1,length(zerocross));
    set( x (x.*last<0) <0) = -1;
    

end