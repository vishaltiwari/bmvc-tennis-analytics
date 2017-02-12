function [graph] = createGraphFromNodeList(node_lst , k_th)
    node_count = size(node_lst,2);
    graph = zeros(node_count,node_count);
    %k_th = 20; % Changed this from 15 -> 20
    for i=1:size(node_lst,2)
        node1 = node_lst(i);
        model1 = node1.model;
        k_min_1 = model1.p1.k;
        k_max_1 = model1.p3.k;
        for j=i+1:size(node_lst,2)
            node2 = node_lst(j);
            if node2.frame - node1.frame > k_th
                break;
            end
%             if j - i > k_th
%                 break;
%             end
            if node1.frame == node2.frame
                continue;
            end
            model2 = node2.model;
            k_min_2 = model2.p1.k;
            k_max_2 = model2.p3.k;

            a1 = model1.compute_a();
            v1_1 = model1.compute_v1(a1);
            a2 = model2.compute_a();
            v1_2 = model2.compute_v1(a2);

            % If they are not overlapping:
            if k_min_2 - k_max_1 > 0
                if k_min_2 - k_max_1 <= k_th
                    min_d = Inf;
                    for k_dash=k_max_1:k_min_2
                        p1 = model1.compute_k_prime(v1_1,a1,k_dash);
                        p2 = model2.compute_k_prime(v1_2,a2,k_dash);
                        d = sqrt((p1(1) - p2(1)) ^2 + (p1(2) - p2(2)) ^2);
                        if d < min_d
                            min_d = d;
                        end
                    end
                    % W between
                    graph(i,j) = min_d;
                end
            else
                % there is overlap:
                flag = 1;
                support_set1 = model1.support_set;
                support_set2 = model2.support_set;
                end_k = min(k_max_1 , k_max_2);
                start_k = max(k_min_1,k_min_2);
                for k_dash=start_k:end_k
                    p1 = support_set1(k_dash - k_min_1 + 1);
                    p2 = support_set2(k_dash - k_min_2 + 1);
                    if ~(p1 == p2)
                        flag = 0;
                        break;
                    end
                    %p1 = model1.compute_k_prime(v1_1,a1,k_dash);
                    %p2 = model2.compute_k_prime(v1_2,a2,k_dash);
                    %d = sqrt((p1(1) - p2(1)) ^2 + (p1(2) - p2(2)) ^2);
                    %if d > 0
                    %    flag = 0;
                    %    break;
                    %end
                end
                if flag == 1
                    graph(i,j) = 10^-100;
                end
            end
        end        
    end
end