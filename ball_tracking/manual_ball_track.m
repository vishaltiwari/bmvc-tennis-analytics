% mark the ball candidates.
addpath('/Neutron9/anurag/AO_2017_segments/ball_track');
seq_no = 20;
video_in = [num2str(seq_no),'_TEST_BALL_DETECTION.avi'];
vidObj = VideoReader(video_in);

vidHeight = vidObj.Height;
vidWidth = vidObj.Width;

mov = struct('cdata',zeros(vidHeight,vidWidth,3,'uint8'),...
    'colormap',[]);

len = 1;
while hasFrame(vidObj)
    mov(len).cdata = readFrame(vidObj);
    len = len+1;
end


ball_points = []; % store its as x, y, frame_nb.
frame_nb = 1;
figure;
while frame_nb < len
    frame = mov(frame_nb).cdata;
    imshow(frame);
    [x, y] = getpts;
    if isempty(x) || isempty(y)
        frame_nb = frame_nb + 1;
        continue;
    end
    ball_points = [ball_points ; x y frame_nb];
    frame_nb = frame_nb + 1;
end

save(['test_data/',num2str(seq_no) , '_ballpoints.mat'] , 'ball_points');