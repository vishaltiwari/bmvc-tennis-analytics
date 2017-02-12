function [length] = getLengthBwExt(start_indx , support_set)
    length = 0;
    for i=start_indx+1:size(support_set,2)
        p1 = support_set(i);
        p2 = support_set(i-1);
        dis = sqrt((p1.x - p2.x)^2 + (p1.y - p2.y)^2);
        length = length + dis;
    end
end