function [] = eventDetection(ballPositions , mov , buffer_len , len , video_out , writeFlag)

    vidOutObj = VideoWriter(video_out);
    open(vidOutObj);

    ball_len = size(ballPositions,2);
    prev_prev_point = ballPositions(1);
    prev_point = ballPositions(2);
    curr_point = ballPositions(3);
    next_point = ballPositions(4);
    next_next_point = ballPositions(5);
    
    prev_a = [(next_point.x + prev_point.x - 2*curr_point.x)/2 , (next_point.y + prev_point.y - 2*curr_point.y)/2];
    prev_v = [(next_point.x-prev_point.x)/2 , (next_point.y-prev_point.y)/2];
    prev_prev_v = [(curr_point.x-prev_prev_point.x)/2 , (curr_point.y-prev_prev_point.y)/2];
    
    events = [];
    
    vertical_v = [prev_v];
    
    for i=4:ball_len
        
        prev_point = curr_point;
        curr_point = next_point;
        next_point = ballPositions(i);
        
        curr_a = [(next_point.x + prev_point.x - 2*curr_point.x)/2 , (next_point.y + prev_point.y - 2*curr_point.y)/2];
        curr_v = [(next_point.x-prev_point.x)/2 , (next_point.y-prev_point.y)/2];
        
        angle_v = acos(dot(curr_v , prev_v) / (norm(curr_v) * norm(prev_v)));

        vertical_v = [vertical_v curr_v(2)];
        % check the y component:
        if curr_v(2) * prev_v(2) < 0
            events = [events curr_point];
        end
        % greater than 15 degrees, then report as an event
%         if angle_v > 90 * (pi/180)
%             disp(abs(norm(curr_v) - norm(prev_v)));
%             if abs(norm(curr_v) - norm(prev_v)) > 5
%                 events = [events curr_point];
%             end
%         end
        
        prev_v = curr_v;
    end
    
    if writeFlag == 1
        event_ptr = 1;
        ball_vel_ptr = 1;
        for i=1+buffer_len: len - buffer_len
            event_point = events(event_ptr);
            if event_point.k == i && event_ptr < size(events,2)
                frame_image = mov(i).cdata;
                if ballPositions(ball_vel_ptr).k + buffer_len == i && ball_vel_ptr < ball_len
                    %frame_image = insertText(frame_image , [ballPositions(ball_vel_ptr).x+5 , ballPositions(ball_vel_ptr).y] , vertical_v(ball_vel_ptr) , 'FontSize' , 18 , 'BoxOpacity',0);
                    ball_vel_ptr = ball_vel_ptr + 1;
                end
                frame_image = insertShape(frame_image , 'circle',[event_point.x event_point.y 5] , 'LineWidth',10);
                writeVideo(vidOutObj , frame_image);
                event_ptr = event_ptr + 1;
            else
                frame_image = mov(i).cdata;
                if ballPositions(ball_vel_ptr).k + buffer_len == i && ball_vel_ptr < ball_len
                    %frame_image = insertText(frame_image , [ballPositions(ball_vel_ptr).x+5 , ballPositions(ball_vel_ptr).y] , vertical_v(ball_vel_ptr) , 'FontSize' , 18 , 'BoxOpacity',0);
                    ball_vel_ptr = ball_vel_ptr + 1;
                end
                writeVideo(vidOutObj , frame_image);
            end
        end
    end
    close(vidOutObj);
end