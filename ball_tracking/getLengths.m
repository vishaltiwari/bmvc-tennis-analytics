function [all_length , node_set] = getLengths(all_paths , node_set)
    all_length = [];
    for i=1:size(all_paths,2)
        path = all_paths(i);
        node_lst=  path{1};
        firstNodeIndx = node_lst(1); %get the first node, and its k_min and k_max
        firstNode = node_set(firstNodeIndx);
        %k_min1 = firstNode.model.p1.k;
        k_max1 = firstNode.model.p3.k;
        len = getLengthBwExt(1,firstNode.model.support_set);
        node_set(firstNodeIndx).model.len = len;
        for j=2:size(node_lst,2)
            nodeIndx = node_lst(j);
            node = node_set(nodeIndx);
            k_min2 = node.model.p1.k;
            k_max2 = node.model.p3.k;
            
            %If there is no-overlap:
            if k_max1 < k_min2
               l = getLengthBwExt(1,node.model.support_set);
               len = len + l;
               node_set(nodeIndx).model.len = l;
               %k_min1 = k_min2;
               k_max1 = k_max2;
               continue;
            end

            %If there is overlap:
            start_k = k_max1;
            l = getLengthBwExt(start_k-k_min2+1,node.model.support_set);
            len = len + l;
            node_set(nodeIndx).model.len = l;
            %k_min1 = k_min2;
            k_max1 = k_max2;
        end
        all_length = [all_length len];
    end
end