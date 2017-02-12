function C = computeModelCost(model,candidate_per_frame,k_start,k_end,d_th)
    a = model.compute_a();
    v1 = model.compute_v1(a);
    
    C = 0;
    for k=k_start:k_end
        pk= model.compute_k_prime(v1,a,k);
        candidates = candidate_per_frame{k};
        for j=1:size(candidates,2)
            if isempty(candidates{j})
                continue;
            end
            % get the distance b/w 
            c = candidates{j};
            dist_sq = (c.x - pk(1))^2 + (c.y-pk(2))^2;
            if dist_sq < d_th^2
                C = C + dist_sq;
            else
                C = C + d_th^2;
            end
        end
    end
end