function [pos] = getBallPositionAtFrame(ballPositions,frame_nb,buffer_len)
    ptr = 1;
    while ptr < size(ballPositions,2)
        %k = ballPositions(ptr).k;
        frame_ball = ballPositions(ptr).k - buffer_len;
        if frame_nb == frame_ball
            pos = ballPositions(ptr);
            return;
        end
        ptr = ptr + 1;
    end
    pos = Point(-1,-1,-1);
end