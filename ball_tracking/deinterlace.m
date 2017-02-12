%video_file = '/media/vishal/0804E66104E6516C/TT_dataset/2012_olympics_men/semi_sample.mp4'; %Table Tennis
%video_file = '/home/vishal/cvit/Table_tennis_dataset/testVideo4.mp4';
function [] = deinterlace(video_file , write_file)
%video_file = '/media/vishal/0804E66104E6516C/Racket_games/Lawn_Tennis/semi_sample.mp4'; % Lawn Tennis

vidObj = VideoReader(video_file);

vidHeight = vidObj.Height;
vidWidth = vidObj.Width;

hdinterlacer = vision.Deinterlacer;

%% Write the de-interlaced data to a file to compare, or to work on later:

%write_file = 'semi_sample_deinterlaced.avi';
v = VideoWriter(write_file);
open(v);


while hasFrame(vidObj)
    image_orig = readFrame(vidObj);
    clearImage = step(hdinterlacer , image_orig);
    writeVideo(v , clearImage);
    %     imshow(clearImage);
    pause(1/vidObj.FrameRate);
end

close(v);
end