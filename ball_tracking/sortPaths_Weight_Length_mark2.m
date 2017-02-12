function [all_path , all_weights , all_length] = sortPaths_Weight_Length_mark2(all_path,all_weights,all_length)
    %Filter the paths which are of length smaller than 20.
    path_threshold = 300;
    all_path = all_path(all_length >= path_threshold);
    all_weights = all_weights(all_length >=path_threshold);
    all_length = all_length(all_length >= path_threshold);
    %Implement the PR algorithm, 
    alpha = 1;
    % Sort all the paths based on the weight and length of each path.
    pathObjectsList = [];
    for i=1:size(all_path,2)
        pathObjectsList = [pathObjectsList PathClass(all_path(i))];
    end
    
end