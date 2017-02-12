%% Video file processing template:
clear;

video_file = 'ov_sample1.avi'; % Lawn Tennis

vidObj = VideoReader(video_file);

vidHeight = vidObj.Height;
vidWidth = vidObj.Width;

while hasFrame(vidObj)
    frame = readFrame(vidObj);
    imshow(frame);
    pause(1/vidObj.FrameRate);
end