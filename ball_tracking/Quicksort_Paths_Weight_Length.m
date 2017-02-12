function [sorted_P, sorted_W, sorted_L] = Quicksort_Paths_Weight_Length(all_path,all_weights,all_length,path_threshold , alpha)
    %Filter the paths which are of length smaller than 20.
    %path_threshold = 200;
    %all_path = all_path(all_length >= path_threshold);
    %all_weights = all_weights(all_length >=path_threshold);
    %all_length = all_length(all_length >= path_threshold);
    %Implement the PR algorithm, 
    %alpha = 1;
    % Sort all the paths based on the weight and length of each path.
    % TODO: implement a q-sort algorithm replacing this bubble sort algo.
    % Use this for comparing two paths, better paths come first.
    % nb = getBetterPath(all_weights(i),all_weights(j),all_length(i),all_length(j),alpha);
    [sorted_P, sorted_W, sorted_L] = qSortPaths(all_path, all_weights , all_length, alpha);

end