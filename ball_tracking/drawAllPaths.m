function [] = drawAllPaths(node_set,mov,allPaths)
    frame_lst = zeros(1,size(allPaths,2));
    for i=1:size(allPaths,2)
        path = allPaths(i);
        node_lst = path{1};
        node = node_set(node_lst(1));
        minFrame = node.model.p1.k;
        frame_lst(i) = minFrame;
    end
    
    % sort the frame_lst, and get the index.
    [B , I] = sort(frame_lst , 'ascend');
    
    %use the I to get the path when iterating paths to show the ball
    %position.
    orderedPaths = [];
    newPathList = [];
    for i=1:size(I,2)
        path = allPaths(I(i));
        newPathList = [newPathList path]; 
        orderedPaths = [orderedPaths path{1}];
    end
    
    findBallTrack([],node_set,mov,orderedPaths);
    
end