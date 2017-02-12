function [all_paths , all_weights] = getAllPairPaths(graph)
    sparse_graph = sparse(graph);
    all_paths = [];
    all_weights = [];
    for i=1:size(graph,2)
        [dist , paths , pred] = graphshortestpath(sparse_graph,i);
        all_paths = [all_paths paths(~cellfun('isempty',paths))];
        all_weights = [all_weights dist(dist~=Inf)];
    end
end