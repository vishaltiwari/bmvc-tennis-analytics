function [] = drawModelsOnImage(node_set,mov)
    %sparse_graph = sparse(graph);
    imshow(mov(10).cdata);
    hold on;
    i=1;
    count = 1;
    while i <= size(node_set,2)
        node_i = node_set(i);
        same_frame_nodes = [node_i];
        for j=i+1:size(node_set,2)
            node_j = node_set(j);
            if node_i.frame == node_j.frame
                same_frame_nodes = [same_frame_nodes node_j];
            else
                break;
            end
        end
        text_flag = 0;
        % plot all the models in same_frame_nodes
        imshow(mov(node_i.frame).cdata);
        hold on;
%         if i >= 2350
%             text_flag = 1;
%             disp('something');
%         end
        for k=1:size(same_frame_nodes,2)
            model = same_frame_nodes(k).model;
            hold on ;
            visualizeTracklet(model,0,0,1);
            if text_flag == 1
                txt1 = num2str(count);
                text(model.p2.x,model.p2.y,txt1)
            end
            count = count  +1;
        end
        i = j;
    end
    
    %[dist , path , pred] = graphshortestpath(graph_sparse,120,1324);
    %[dist , paths , pred] = graphshortestpath(sparse_graph,1);
    
    %path contains the nodes, which are the tracklets. Visualize them as
    %such
%     for i=1:size(paths,2)
%         model = node_set(path(i));
%         hold on;
%         visualizeTracklet(model,0,0,1); %use the bounds of the trackets
%     end
end