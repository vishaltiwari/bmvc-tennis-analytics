function [ballPositions] = writeBestPath(node_set , mov , len, buffer_len ,bestP, video_out)
    
    ballPositions = [];
    vidOutObj = VideoWriter(video_out);
    open(vidOutObj);
    
    % write out to a video to trace the ball path:
    court_image = mov(buffer_len+1).cdata;
    videOutObjImage = VideoWriter('ball_trace.avi');
    open(videOutObjImage);

    path = bestP{1};
    
    total_nodes = size(path,2);
    
    node_ptr = 1;
    
    start_k = node_set(path(node_ptr)).model.p1.k + buffer_len;
    end_k = node_set(path(node_ptr)).model.p3.k + buffer_len;
    curr_model = node_set(path(node_ptr)).model;
    
    start_k_next = node_set(path(node_ptr+1)).model.p1.k + buffer_len;
    next_model = node_set(path(node_ptr+1)).model;
    
    curr_support_set = node_set(path(node_ptr)).model.support_set;
    
    %loop over the video, i is the frame nb.
    i = 1+buffer_len;
    while i <= len-buffer_len
        if i<start_k
            writeVideo(vidOutObj , mov(i).cdata);
            i = i + 1;
            continue
        end
        % plot the current model position.
        if i <= end_k
            ball_pos = curr_support_set(i - start_k + 1);
        else
            if i >= start_k_next && node_ptr+1 <= total_nodes
                node_ptr = node_ptr + 1;
                
                start_k = node_set(path(node_ptr)).model.p1.k + buffer_len;
                end_k = node_set(path(node_ptr)).model.p3.k + buffer_len;
                curr_model = node_set(path(node_ptr)).model;
                
                curr_support_set = node_set(path(node_ptr)).model.support_set;
                if node_ptr+1 <= total_nodes
                    start_k_next = node_set(path(node_ptr+1)).model.p1.k + buffer_len;
                    end_k_next = node_set(path(node_ptr+1)).model.p1.k + buffer_len;
                    next_model = node_set(path(node_ptr+1)).model;
                end
                %end_k_next = node_set(path(node_ptr+1)).model.p3.k + buffer_len;
                
                if i > end_k && (node_ptr + 1) <= size(path,2);
                    node_ptr = node_ptr + 1;
                    start_k = node_set(path(node_ptr)).model.p1.k + buffer_len;
                    end_k = node_set(path(node_ptr)).model.p3.k + buffer_len;
                    curr_model = node_set(path(node_ptr)).model;
                    
                    curr_support_set = node_set(path(node_ptr)).model.support_set;
                    if node_ptr+1 <= total_nodes
                        start_k_next = node_set(path(node_ptr+1)).model.p1.k + buffer_len;
                        end_k_next = node_set(path(node_ptr+1)).model.p1.k + buffer_len;
                        next_model = node_set(path(node_ptr+1)).model;
                    end
                    continue;
                else
                    if (i - start_k + 1) <= size(curr_support_set,2)
                        ball_pos = curr_support_set(i - start_k + 1);
                    end
                end
                    
            elseif i < start_k_next
                % we will have to interpolate.
                prev_a = curr_model.compute_a();
                prev_v = curr_model.compute_v1(prev_a);
                
                next_a = next_model.compute_a();
                next_v = next_model.compute_v1(next_a);
                
                min_d = intmax; %min value will be the point of intersection.
                min_k = -1;
                for k=end_k:start_k_next
                    pos_prev = curr_model.compute_k_prime(prev_v,prev_a,k-buffer_len);
                    pos_next = next_model.compute_k_prime(next_v,next_a,k-buffer_len);
                    d = sqrt((pos_prev(1) - pos_next(1))^2 + (pos_prev(2) - pos_next(2))^2);
                    if d < min_d
                        min_d = d;
                        min_k = k;
                    end
                end
                
                while i <= min_k
                    pos_prev = curr_model.compute_k_prime(prev_v,prev_a,i-buffer_len);
                    frame_image = mov(i).cdata;
                    %frame_image = insertShape(frame_image , 'circle' , [pos_prev(1),pos_prev(2),5] , 'color' , 'yellow');
                    frame_image = insertShape(frame_image , 'rectangle' , [pos_prev(1)-5,pos_prev(2)-5,10,10] , 'color' , 'yellow');
                    court_image = insertShape(court_image, 'circle' , [pos_prev(1),pos_prev(2),2] , 'color' , 'black');
                    writeVideo(videOutObjImage , court_image);
                    writeVideo(vidOutObj , frame_image);
                    ball_point = Point(pos_next(1) , pos_next(2) , i-buffer_len);
                    ballPositions = [ballPositions ball_point];
                    i = i +1;
                end
                
                while i < start_k_next
                    pos_next = next_model.compute_k_prime(next_v,next_a,i-buffer_len);
                    frame_image = mov(i).cdata;
                    %frame_image = insertShape(frame_image , 'circle' , [pos_next(1),pos_next(2),5] , 'color' , 'yellow');
                    frame_image = insertShape(frame_image , 'rectangle' , [pos_next(1)-5,pos_next(2)-5,10,10] , 'color' , 'yellow');
                    court_image = insertShape(court_image, 'circle' , [pos_next(1),pos_next(2),2] , 'color' , 'black');
                    writeVideo(videOutObjImage , court_image);
                    writeVideo(vidOutObj , frame_image);
                    ball_point = Point(pos_next(1) , pos_next(2) , i-buffer_len);
                    ballPositions = [ballPositions ball_point];
                    i = i +1;
                end
                continue;
                %writeVideo(vidOutObj , mov(i).cdata);
                %ball_pos = Point(-1,-1,-1);
            end
        end
        
        frame_image = mov(i).cdata;
        %frame_image = insertShape(frame_image , 'circle' , [ball_pos.x,ball_pos.y,5] , 'color' , 'red');
        frame_image = insertShape(frame_image , 'rectangle' , [ball_pos.x-5,ball_pos.y-5,10,10] , 'color' , 'red');
        court_image = insertShape(court_image, 'circle' , [ball_pos.x,ball_pos.y,2] , 'color' , 'yellow');
        writeVideo(videOutObjImage , court_image);
        ball_point = Point(ball_pos.x , ball_pos.y , i-buffer_len);
        ballPositions = [ballPositions ball_point];
        writeVideo(vidOutObj , frame_image);
        i = i +1;
    end
    close(vidOutObj);
    close(videOutObjImage);
end