function graph = createGraphFromList(models_lst)
    model_count = size(models_lst,2);
    graph = zeros(model_count,model_count);
    k_th = 15;
    N = 7;
    for i=1:size(models_lst,2)
        model1 = models_lst(i);
        k_min_1 = model1.p1.k;
        k_max_1 = model1.p3.k;
        for j=i+1:size(models_lst,2)
            model2 = models_lst(j);
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
                    graph(i,j) = min_d;
                else
                    graph(i,j) = Inf;
                end
            else
                % there is overlap:
                flag = 1;
                for k_dash=k_min_2:k_max_1
                    d1 = model1.compute_k_prime(v1_1,a1,k_dash);
                    d2 = model2.compute_k_prime(v1_2,a2,k_dash);
                    d = sqrt((d1(1) - d2(1)) ^2 + (d1(2) - d2(2)) ^2);
                    if d > 0
                        flag = 0;
                        break;
                    end
                end
                if flag == 1
                    graph(i,j) = 0;
                else
                    graph(i,j) = Inf;
                end
            end
        end
    end
end