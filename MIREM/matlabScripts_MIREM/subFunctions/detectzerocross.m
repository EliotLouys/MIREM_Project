function [zerocross, set]  = detectzerocross(x)
    idx                    = [1:length(x)];
    last                   = circshift(x,1);
    zerocross              = idx(x.*last<0);  %  = TRUE when x crosses zero
    set                    = ones(1,length(zerocross));
    set( x (x.*last<0) <0) = -1;
    

end