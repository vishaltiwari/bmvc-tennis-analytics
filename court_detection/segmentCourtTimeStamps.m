function [] = segmentCourtTimeStamps(video_in , video_write , template_in)
%Need to load/setup vl_feat before running this.

%video_in = '/home/vishal/cvit/tennis_info_sys/ball_tracking/sample_data/australian_open_sample1.mp4';
vidObj = VideoReader(video_in,'Motion JPEG AVI');

%video_write = 'detect_Lawn_Tennis_using_template_1_video2.avi';
v = VideoWriter(video_write);
%v.VideoCompressionMethod = vidObj.VideoCompressionMethod;
v.FrameRate = vidObj.FrameRate;


open(v);

%Get the features of the template provided.
[templateFeatures , templateDescriptors , court_pos] = getTemplateFeaturesDescriptors_mark2(template_in);

while hasFrame(vidObj)
    curr_image = readFrame(vidObj);
    
    %Get the features of the current frame.
    [frameFeatures , frameDescriptors] = getFrameFeatureDescriptor_mark2(curr_image);
    
    %Match template desciptors to the current frame descriptoes.
    [court_flag , tform] = matchDescriptors(templateFeatures , templateDescriptors , frameFeatures , frameDescriptors,curr_image , template_in);
    if court_flag == 0
        tempImage = insertText(curr_image,[5 5],'No court','BoxColor','red','BoxOpacity',0.4,'TextColor','white');
        %continue
    else
        posNew = transformPointsForward(tform , court_pos);
        p = [posNew 3*ones(size(posNew,1),1)];
        tempImage = insertShape(curr_image,'circle',p , 'LineWidth',5);
        tempImage = insertText(tempImage,[5 5],'court found','BoxColor','green','BoxOpacity',0.4,'TextColor','white');
    end
    %figure ; imshow(tempImage);
    writeVideo(v,tempImage);
%    pause(1/vidObj.FrameRate);
end
close(v);
%use the flag_lst to segment out the video into play and non-play.

end