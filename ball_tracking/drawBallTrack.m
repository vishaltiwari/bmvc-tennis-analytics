function [] = drawBallTrack(graph,node_set,mov,best_paths)

    sparse_graph = sparse(graph);
    imshow(mov(10).cdata);
    hold on;
    %[dist1 , path1 , pred1] = graphshortestpath(sparse_graph,1476,2099); % Path1
    
    %concat a list:
    
    %[dist , paths , pred] = graphshortestpath(sparse_graph,1);
%    path = [path1 path2 path3 path4 path5 path6 path7];
    %path = [path1];
    %path contains the nodes, which are the tracklets. Visualize them as
    %such
    for j=1:size(best_paths,2)
        p = best_paths(j);
        path = p{1};
        %imshow(mov(10).cdata);
        %hold on;
        for i=1:size(path,2)
            model = node_set(path(i)).model;
            hold on;
            visualizeTracklet(model,0,0,1); %use the bounds of the trackets
        end
    end
end