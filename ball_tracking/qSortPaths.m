function [paths, weights, lengths] = qSortPaths(paths, weights, lengths,alpha)
    
    i = 1;
    j = size(weights,2);
    if i==j
        paths = paths(i);
        weights = weights(i);
        lengths = lengths(i);
        return;
    end

    pivot= j;
    mid = floor((i + j) /2);
    
    % Swappp mid, and last path
    tmp_p = paths(pivot);
    tmp_w = weights(pivot);
    tmp_l = lengths(pivot);
    paths(pivot) = paths(mid);
    weights(pivot) = weights(mid);
    lengths(pivot) = lengths(mid);
    paths(mid) = tmp_p;
    weights(mid) = tmp_w;
    lengths(mid) = tmp_l;
    j = j-1;
    while i <= j
        while i < j
            nb = getBetterPath(weights(pivot),weights(i),lengths(pivot),lengths(i),alpha);
            if nb == 2 || nb==0
                i = i+1;
                continue;
            else
                break;
            end
        end
        
        while i < j
            nb = getBetterPath(weights(pivot),weights(j),lengths(pivot),lengths(j),alpha);
            if nb == 1 || nb==0
                j = j-1;
                continue;
            else
                break;
            end
        end
        
        %Swap paths at i,j.
        tmp_p = paths(i);
        tmp_w = weights(i);
        tmp_l = lengths(i);
        paths(i) = paths(j);
        weights(i) = weights(j);
        lengths(i) = lengths(j);
        paths(j) = tmp_p;
        weights(j) = tmp_w;
        lengths(j) = tmp_l;
        
        if i>=j
            break;
        end
        
        i = i+1;
        j = j-1;
        
    end
    
    %move the pivot back to j.
    tmp_p = paths(j);
    tmp_w = weights(j);
    tmp_l = lengths(j);
    paths(j) = paths(pivot);
    weights(j) = weights(pivot);
    lengths(j) = lengths(pivot);
    paths(pivot) = tmp_p;
    weights(pivot) = tmp_w;
    lengths(pivot) = tmp_l;
    
    [path_left, weights_left , length_left] = qSortPaths(paths(1:j) , weights(1:j) , lengths(1:j) , alpha);
    [path_right, weights_right , length_right] = qSortPaths(paths(j+1:pivot), weights(j+1:pivot), lengths(j+1:pivot) , alpha);
    
    paths = [path_left path_right];
    weights = [weights_left weights_right];
    lengths = [length_left length_right];
    
end