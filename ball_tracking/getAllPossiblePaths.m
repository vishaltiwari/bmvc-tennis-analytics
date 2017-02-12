function [paths] = getAllPossiblePaths(graph)
    addpath('/home/vishal/libraries/pathbetweennodes-pkg/pathbetweennodes');
    for i=1:size(graph,2)
        for j=i+1:size(graph,2)
            path = pathbetweennodes(graph , i , j,true);
        end
    end
end