function [] = createGraph(tracklets_arr,model_count)
    graph = zeors(model_count,model_count);
    k_th = 15;
    N = 7;
    model_num=0;
    for i=1:size(tracklets_arr,2)
        model_lst = tracklets_arr{i}; %get the track models at frame i;
        for j=1:size(model_lst,2)
            model1 = model_lst(i);
            model_num = model_num + 1;
            model1_index = model_num;
            k_min_1 = model1.p1.k;
            k_max_1 = model1.p3.k;
            for l=1:2*N+1
                ahead_ptr= i+l;
                ahead_model_lst = tracklets_arr{ahead_ptr};
                for q=1:size(ahead_model_lst,2)
                    model2 = ahead_model_lst(q);
                    model_num = model_num + 1;
                    model2_index = model_num;
                    k_min_2 = model2.p1.k;
                    k_max_2 = model2.p3.k;
                    
                    a1 = model1.compute_a();
                    v1_1 = model1.compute_v1(a1);
                    a2 = model2.compute_a();
                    v1_2 = model2.compute_v1(a2);
                    
                    % If they are not overlapping:
                    if k_min_2 - k_max_1 > 0
                        if k_min_2 - k_max_1 < k_th
                            min_d = 1000; % NOTE: HARD CODED value. should be max(width,height) of image.
                            for k_dash=k_max_1:k_min_2
                                d1 = model1.compute_k_prime(v1_1,a1,k_dash);
                                d2 = model2.compute_k_prime(v1_2,a2,k_dash);
                                d = sqrt((d1(1) - d2(1)) ^2 + (d1(2) - d2(2)) ^2);
                                if d < min_d
                                    min_d = d;
                                end
                            end
                            % W between 
                            graph(model1_index,model2_index) = min_d;
                        else
                            graph(model1_index,model2_index) = Inf;
                        end
                    else
                        % If there is overlap:
                        flag = 1;
                        for k_dash=k_min_2:k_max_1
                            d1 = model1.compute_k_prime(v1_1,a1,k_dash);
                            d2 = model2.compute_k_prime(v1_2,a2,k_dash);
                            d = sqrt((d1(1) - d2(1)) ^2 + (d1(2) - d2(2)) ^2);
                            if d >= 1
                                flag = 0;
                                break;
                            end
                        end
                        if flag == 1
                            graph(model1_index,model2_index) = 0;
                        else
                            graph(model1_index,model2_index) = Inf;
                        end
                    end
                end
            end
        end
    end
end