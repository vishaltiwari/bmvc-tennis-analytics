function [all_length]= getAllLengthPaths(all_paths , node_set)
    all_length = [];
    for i=1:size(all_paths,2)
        path = all_paths(i);
        
        node_lst = path{1};
        
        n = node_set(1);
        min_k = n.model.p1.k;
        max_k = n.model.p3.k;
        l=0;
        for node_ptr=1:size(node_lst,2)
            node = node_set(node_ptr);
            k_min = node.model.p1.k;
            k_max = node.model.p3.k;

            %If there is no overlap, we need to reset min_k,max_k
            if min_k > max_k
                l = l + max_k - min_k + 1;
                min_k = Inf;
                max_k = -1*Inf;
            end
            
            if k_min < min_k
                min_k = k_min;
            end
            if k_max > max_k
                max_k = k_max;
            end
        end
        
        l = l + max_k - min_k + 1;
        all_length = [all_length l];

    end
end