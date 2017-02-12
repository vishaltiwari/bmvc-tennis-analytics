function [flag_compatibility] = isTrackCompatible(path1,path2,node_set)
    i=1;
    j=1;
    flag_compatibility = 1;
    node_lst1 = path1{1};
    node_lst2 = path2{1};
    while i <= size(node_lst1,2) && j<=size(node_lst2,2)
        nodei_index = node_lst1(i);
        nodei = node_set(nodei_index);
        k_min1 = nodei.model.p1.k;
        k_max1 = nodei.model.p2.k;
        %support_set1 = nodei.model.support_set;
        nodej_index = node_lst2(j);
        nodej = node_set(nodej_index);
        k_min2 = nodej.model.p1.k;
        k_max2 = nodej.model.p3.k;
        %case when there is no intersection:
        if k_max1 < k_min2
            %increment i;
            i = i+1;
            continue;
        end
        if k_min1 > k_max2
            %inc j
            j = j+1;
            continue;
        end

        %Else there is an intersection:
        support_set1 = nodei.model.support_set;
        support_set2 = nodej.model.support_set;
        ss1_size = size(support_set1,2);
        ss2_size = size(support_set2,2);
        p=1;
        q=1;
        while support_set1(p).k < support_set2(q).k && p <= ss1_size
            p = p + 1;
        end

        while support_set2(q).k < support_set1(p).k && q <= ss2_size
            q = q + 1;
        end

        % now the k for both the support set are the same.
        while p <= ss1_size && q <= ss2_size
            point1 = support_set1(p);
            point2 = support_set2(q);
            if point1 == point2
                flag_compatibility = 0;
                break;
            end
            p = p + 1;
            q = q + 1;
        end
        if flag_compatibility ==0
            break;
        end
        if p >= ss1_size && q >= ss2_size
            i = i+1;
            j = j+1;
        elseif p >= ss1_size
            i = i+1;
        elseif q >= ss2_size
            j = j+1;
        end
        
    end
end