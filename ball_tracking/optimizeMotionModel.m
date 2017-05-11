function optimizedModelsList = optimizeMotionModel(model_lst,candidate_per_frame,k_start,k_end , m_th , d_th , frame_img, vidOutObj , buffer_len, writeFlag)
    optimizedModelsList = [];
    %m_th = 4; % Added a new filtering parameter.
    %d_th = 5; % defined in the paper, changed from 3 -> 5
    dth_sq = d_th ^ 2;
    for i=1:size(model_lst,2)
        
        model = model_lst(i);
        if writeFlag == 1
            frame_img = visualizeTracklet(frame_img, model,candidate_per_frame,k_start,k_end);
            writeVideo(vidOutObj , frame_img);
        end
        % move in the frame direction, k_min and k_max;
        % This starts of with the seed pixels
        k_min = model.p1.k;
        k_max = model.p3.k;
        C_prev = computeModelCost(model,candidate_per_frame,k_start,k_end,d_th);
        support_set = [model.p1 model.p2 model.p3];
        while true
            a = model.compute_a();
            v1 = model.compute_v1(a);
            %move back in frames:
            k = k_min-1;
            found_flag = 1;
            tmp_k_min_point = model.p1;
            % Going back in frame.
            while k > k_start
                % for the specific k
                if found_flag == 0
                    break;
                end
                candidates = candidate_per_frame{k};
                predicted_pos = model.compute_k_prime(v1,a,k);
                found_flag = 0;
                for c=1:size(candidates,2)
                    if isempty(candidates{c})
                        continue;
                    end
                    % get the distance between predicted and all
                    % candidates:
                    dis = (predicted_pos(1) - candidates{c}.x)^2 +  (predicted_pos(2) - candidates{c}.y)^2;
                    if dis < dth_sq
                        %found a new support point.
                        tmp_k_min_point = Point(candidates{c}.x , candidates{c}.y , candidates{c}.k);
                        found_flag = 1;
                        support_set = [tmp_k_min_point support_set];
                    end
                end
                k = k - 1;
            end
            
            tmp_k_max_point = model.p3;
            % go
            k = k_max+1;
            found_flag=1;
            while k < k_end
                % for the specific k
                if found_flag == 0
                    break;
                end
                candidates = candidate_per_frame{k};
                predicted_pos = model.compute_k_prime(v1,a,k);
                found_flag = 0;
                for c=1:size(candidates,2)
                    if isempty(candidates{c})
                        continue;
                    end
                    % get the distance between predicted and all
                    % candidates:
                    dis = (predicted_pos(1) - candidates{c}.x)^2 +  (predicted_pos(2) - candidates{c}.y)^2;
                    if dis < dth_sq
                        %found a new support point.
                        tmp_k_max_point = Point(candidates{c}.x , candidates{c}.y , candidates{c}.k);
                        found_flag = 1;
                        support_set = [support_set tmp_k_max_point];
                    end
                end
                k = k + 1;
            end
            
            % End the loop, when there is no change in k_min, and K-max
            % points.
            if tmp_k_max_point == model.p3 && tmp_k_min_point == model.p1
                break;
            end
            
            % If reached here, means we need to update the model, according
            % to the new k_min and K_max, and find the k_mid.
            % find the mid point:
            k_max = tmp_k_max_point.k;
            k_min = tmp_k_min_point.k;
            k_mid = 0;
            min_value = k_end - k_start + 1;
            for k=k_min+1:k_max-1
                predicted_pos = model.compute_k_prime(v1,a,k);
                candidates = candidate_per_frame{k};
                for c=1:size(candidates,2)
                    if isempty(candidates{c})
                        continue;
                    end
                    dis_sq = (predicted_pos(1) - candidates{c}.x)^2 +  (predicted_pos(2) - candidates{c}.y)^2;
                    if dis_sq < dth_sq
                        val = abs((k_max - k) - (k-k_min));
                        if val <= min_value
                            min_value = val;
                            k_mid = candidates{c};
                        end
                    end
                end
            end

            %Create a new model, and find the cost of the model:
            model_new = Model(tmp_k_min_point,k_mid,tmp_k_max_point);
            model_new.support_set = support_set;
            % calculate the cost, and if > prev_cost, break;
            C = computeModelCost(model_new,candidate_per_frame,k_start,k_end,d_th);
            if C > C_prev
                break;
            end
            C_prev = C;

            % Update the model, with p1= k_min, p2=k_mid, p3 = k_max
            %model.updateModel(tmp_k_min_point , k_mid , tmp_k_max_point);
            model = model_new;
            if writeFlag == 1
                frame_img = visualizeTracklet(frame_img, model,candidate_per_frame, k_start,k_end);
                writeVideo(vidOutObj , frame_img);
            end
            %visualizeTracklet(model,k_start,k_end,1);
            %scatter([model.p1.x ; model.p2.x ; model.p3.x] , [model.p1.y ; model.p2.y ; model.p3.y],100,'d');
        end
        %If the models has less than m_th support_set, don't add it.
        if size(model.support_set,2) >= m_th
            optimizedModelsList = [optimizedModelsList model];
        end
    end
end