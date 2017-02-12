function [best_path_set , best_path_weight, best_path_length] = getFinalTracks_APSP(all_path,all_weights,all_length,node_set)
    
    % Iterate over the sorted path
    best_path_set = [all_path(1)];
    best_path_length = [all_length(1)];
    best_path_weight = [all_weights(1)];
    for i=2:size(all_path,2)
        path = all_path(i);
        flag = 1;
        for j=1:size(best_path_set,2)
            path_B = best_path_set(j);
            flag_comp = isTrackCompatible(path , path_B,node_set);
            %flag_comp = isTrackCompatible_mark2(path , path_B,node_set); % Just testing this out.
            if flag_comp == 0
                flag = 0;
                break;
            end
        end
        if flag == 1
            best_path_set = [best_path_set path];
            best_path_length = [best_path_length all_length(i)];
            best_path_weight = [best_path_weight all_weights(i)];
        end
    end
end