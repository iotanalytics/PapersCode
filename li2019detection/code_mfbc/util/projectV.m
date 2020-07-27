function [T] = projectV(T, V)


if length(V) > 2

    if ~issorted(V)
        error('V must be sorted in ascending order')
    end
    
    edges = [];
    
    for i=1:(length(V)-1)
        edges = [edges  0.5 * (V(i) + V(i+1))];
    end
    
    edges = [-inf edges inf];
    [~, bin] = histc(T(:), edges);
    %T(T > 1) = 1;
    %T(T < 0) = 0;
    
    T = reshape(V(bin), size(T));
    
else
    
    T = (V(2) - V(1)) * double(T > (mean(V))) + V(1);
    
end


end

