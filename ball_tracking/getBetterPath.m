function [path] = getBetterPath(w1,w2,l1,l2,alpha)
    path = 1;
    if w1-w2 < alpha * (l1-l2)
        path = 1;   %Path 1 is better
    elseif w1-w2 > alpha * (l1-l2)
        path = 2;   %Path 2 is better
    else
        if w1 < w2
            path = 1;
        elseif w1 > w2
            path = 2;
        else
            if l1 < l2
                path = 2;
            elseif l1 > l2
                path = 1;
            else
                path = 0;
            end
        end
        path = 0;   %Both paths are equal
    end
end