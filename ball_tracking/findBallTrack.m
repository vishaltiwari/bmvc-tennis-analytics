function [] = findBallTrack(graph,node_set,mov,bestP)

    video_out = './TRGMCOutputVideo/ov_aus_open_sample1_BALL_PATH_APSP.avi';

    vidOutObj = VideoWriter(video_out);
    open(vidOutObj);
    
    sparse_graph = sparse(graph);
    imshow(mov(10).cdata);
    hold on;
%     [dist1 , path1 , pred1] = graphshortestpath(sparse_graph,135,259); % Path1
%     [dist2 , path2 , pred2] = graphshortestpath(sparse_graph,279,616); % Path2
%     [dist3 , path3 , pred3] = graphshortestpath(sparse_graph,738,1073); % Path3
%     [dist4 , path4 , pred4] = graphshortestpath(sparse_graph,1152,1607); % Path4
%     [dist5 , path5 , pred5] = graphshortestpath(sparse_graph,1674,2033); % Path5
%     [dist6 , path6 , pred6] = graphshortestpath(sparse_graph,2061,2350); % Path6
%     [dist7 , path7 , pred7] = graphshortestpath(sparse_graph,2413,2480); % Path6
    %[dist , paths , pred] = graphshortestpath(sparse_graph,1);
    %path = [path1 path2 path3 path4 path5 path6 path7];
    %path = [path1];
    %path contains the nodes, which are the tracklets. Visualize them as
    %such
%     for i=1:size(path,2)
%         model = node_set(path(i)).model;
%         hold on;
%         visualizeTracklet(model,0,0,1); %use the bounds of the trackets
%     end
    i=1;
    model_curr = node_set(bestP(i)).model;
    model_next = node_set(bestP(i+1)).model;
    k_min_curr = model_curr.p1.k;
    k_max_curr = model_curr.p3.k;
    ss_ptr = 1;
    curr_support_set = model_curr.support_set;
    
    k_min_next = model_next.p1.k;
    k_max_next = model_next.p3.k;
    next_support_set = model_next.support_set;
    
    a = model_curr.compute_a();
    v1 = model_curr.compute_v1(a);
    for k=7:size(mov,2)
        curr_img = mov(k+6).cdata;
        if k < k_min_curr
            %not found, leave it, continue
            writeVideo(vidOutObj,curr_img);
            imshow(curr_img);
            hold on ;
            continue;
        end
        
        % If there is overlap with the next model, change the model to the
        % next model.
        if k==k_min_next
            i = i+1;
            model_curr = model_next;
            k_min_curr = model_curr.p1.k;
            k_max_curr = model_curr.p3.k;
            a = model_curr.compute_a();
            v1 = model_curr.compute_v1(a);
            
            model_next = node_set(bestP(i+1)).model;
            ss_ptr = 1;
        end
        
        %p_dash = model_curr.compute_k_prime(v1,a,k);
        p_dash = curr_support_set(ss_ptr);
        ss_ptr = ss_ptr + 1;
        curr_img = insertShape(curr_img, 'circle' ,[p_dash(1) , p_dash(2) , 5] );
        imshow(curr_img);
        writeVideo(vidOutObj,curr_img);
        hold on ;

        
        %If model lasts till the end of the model completion.
%         if k == k_max_curr
%             % Get the next path:
%             i = i+1;
%             model_curr = node_set(bestP(i)).model;
%             k_min_curr = model_curr.p1.k;
%             k_max_curr = model_curr.p3.k;
%             a = model_curr.compute_a();
%             v1 = model_curr.compute_v1(a);
%         end
    end
    close(vidOutObj);
end