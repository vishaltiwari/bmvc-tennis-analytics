function [flag_compatibility] = isTrackCompatible_mark2(path1,path2,node_set)
    % Just use the overlap of the time frames, if not compatible.
    node1_lst = path1{1};
    node1 = node1_lst(1);
    node1Last = node1_lst(size(node1_lst,2));
    k_min1 = node_set(node1).model.p1.k;
    k_max1 = node_set(node1Last).model.p3.k;
    
    node2_lst = path2{1};
    node2 = node2_lst(1);
    node2Last = node2_lst(size(node2_lst,2));
    k_min2 = node_set(node2).model.p1.k;
    k_max2 = node_set(node2Last).model.p3.k;
    
    flag_compatibility = 0;
    
    if k_min2 >= k_max1 || k_min1 >= k_max2
        flag_compatibility = 1;
    end
    
end