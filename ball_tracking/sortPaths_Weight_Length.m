function [all_path , all_weights , all_length] = sortPaths_Weight_Length(all_path,all_weights,all_length)
    %Filter the paths which are of length smaller than 20.
    path_threshold = 300;
    all_path = all_path(all_length >= path_threshold);
    all_weights = all_weights(all_length >=path_threshold);
    all_length = all_length(all_length >= path_threshold);
    %Implement the PR algorithm, 
    alpha = 1;
    % Sort all the paths based on the weight and length of each path.
    %TODO: implement a q-sort algorithm replacing this bubble sort algo.
    for i=1:size(all_path,2)
        for j=i+1:size(all_path,2)
            nb = getBetterPath(all_weights(i),all_weights(j),all_length(i),all_length(j),alpha);
            if nb == 2
                %Swap
                tmp_p = all_path(i);
                tmp_w = all_weights(i);
                tmp_l = all_length(i);
                all_path(i) = all_path(j);
                all_weights(i) = all_weights(j);
                all_length(i) = all_length(j);
                all_path(j) = tmp_p;
                all_weights(j) = tmp_w;
                all_length(j) = tmp_l;
            end
        end
    end
end