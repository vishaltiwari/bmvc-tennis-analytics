function [] = drawSingleBallTrack(graph,node_set,mov,path1)

    sparse_graph = sparse(graph);
    imshow(mov(10).cdata);
    hold on;
    %[dist1 , path1 , pred1] = graphshortestpath(sparse_graph,135,259); % Path1 % Path1
    % 1159 - 1652
    
    for i=1:size(path1,2)
        model = node_set(path1(i)).model;
        hold on;
        visualizeTracklet(model,0,0,1); %use the bounds of the trackets
    end
    
end