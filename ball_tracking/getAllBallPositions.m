function [ballPositions] = getAllBallPositions(node_set,mov,allPaths,video_out)

    %video_out = './TRGMCOutputVideo/ov_aus_open_sample1_BALL_PATH_APSP.avi';

    vidOutObj = VideoWriter(video_out);
    open(vidOutObj);

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
    model_lst = [];
    for i=1:size(I,2)
        path = allPaths(I(i));
        model_lst = [model_lst path{1}];
    end
    
    % f1 => 3 - 12.
    % f2 => 12 -> 17.
    %Iterating of the frame numbers, and get the best model from that
    %frame.
    i=1;
    ballPositions = [];%zeros(1,size(mov,2));
    buffer_len = 6; % Same as the one used padding ones on each size of the mov object.
    bestNodeIndx = model_lst(i);
    bestNode = node_set(bestNodeIndx);
    l1 = bestNode.model.len;
    f = 1;
    for f=1:size(mov,2)
        %look for the models in a buffer range of 7.
        if f+buffer_len > size(mov,2)
            break;
        end
        for j=1:size(model_lst,2)
            node = node_set(model_lst(j));
            kmin = node.model.p1.k;
            kmax = node.model.p3.k;
            if f < kmin || f > kmax
                continue;
            end
            
            %If current f is not in the range of the model.
            if f > bestNode.model.p3.k
                bestNodeIndx = model_lst(j);
                bestNode = node_set(bestNodeIndx);
                l1 = bestNode.model.len;
                i=j;
                continue;
            end
            
            if f >= node.model.p1.k && f <= node.model.p3.k
                l2 = node.model.len;
                if l2 > l1
                    bestNodeIndx = model_lst(j);
                    bestNode = node_set(bestNodeIndx);
                    l1 = bestNode.model.len;
                    i=j;
                end
            end
        end
        %From the bestNode, get the 
        support_set = bestNode.model.support_set;
        k_min = bestNode.model.p1.k;
        if f < k_min
            %imshow(mov(f + buffer_len).cdata);
            writeVideo(vidOutObj,mov(f + buffer_len).cdata);
            %hold on;
            continue;
        end
        pred_flag = 0;
        if (f-k_min+1) > size(support_set,2)
            pred_flag = 1;
            %imshow(mov(f + buffer_len).cdata);
            writeVideo(vidOutObj,mov(f + buffer_len).cdata);
            %hold on;
            %continue;
            model_curr = bestNode.model;
            a = model_curr.compute_a();
            v1 = model_curr.compute_v1(a);
            p = model_curr.compute_k_prime(v1,a,f);
            point = Point(p(1),p(2),f);
        else
           point = support_set(f-k_min+1); 
        end
        ballPositions = [ballPositions point]; % location of the ball at diff frames.
        if pred_flag == 1
            curr_img = insertShape(mov(f + buffer_len).cdata, 'circle' ,[point.x , point.y , 5], 'Color','yellow');
        else
            curr_img = insertShape(mov(f + buffer_len).cdata, 'circle' ,[point.x , point.y , 5], 'Color','red');
        end
        
        writeVideo(vidOutObj,curr_img);
        %imshow(curr_img);
        %hold on;
    end
    close(vidOutObj);
end