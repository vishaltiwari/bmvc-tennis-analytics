function [] = videoGlobalMotionCompensation(video_in , video_write , template_in)
%Need to load/setup vl_feat before running this.
tic
%video_in = '/home/vishal/cvit/tennis_info_sys/ball_tracking/sample_data/australian_open_sample1.mp4';
vidObj = VideoReader(video_in);

%video_write = 'detect_Lawn_Tennis_using_template_1_video2.avi';
%video_write_two = strcat(video_write , 'court_detection');
v = VideoWriter(video_write,'Motion JPEG AVI');
%v2 = VideoWriter(video_write_two);
%v.Quality = 40;
%v2.Quality = 40;
open(v);
%open(v2);

%Get the features of the template provided.
[templateFeatures , templateDescriptors , court_pos] = getTemplateFeaturesDescriptors_mark2(template_in);
base_image = imread(template_in);
%height = vidObj.Height;
%width = vidObj.Width;
hdinterlacer = vision.Deinterlacer;
while hasFrame(vidObj)
    curr_image = readFrame(vidObj);

    %Get the features of the current frame.
    [frameFeatures , frameDescriptors] = getFrameFeatureDescriptor_mark2(curr_image);
    
    %Match template desciptors to the current frame descriptoes.
    [court_flag , tform] = matchDescriptors(templateFeatures , templateDescriptors , frameFeatures , frameDescriptors,curr_image , template_in);
    if court_flag == 0
        %tempImage = insertText(curr_image,[15 15],'No court','BoxColor','red','BoxOpacity',0.4,'TextColor','white');
        continue
    else
        invtform = invert(tform);
        % deinterlace the frame
        curr_image = step(hdinterlacer , curr_image);
    
        image_new = imwarp(curr_image , invtform , 'OutputView',imref2d(size(base_image)));
        writeVideo(v,image_new);
        
        %posNew = transformPointsForward(tform , court_pos);
        %p = [posNew 3*ones(size(posNew,1),1)];
        %tempImage = insertShape(curr_image,'circle',p , 'LineWidth',5);
        %tempImage = insertText(tempImage,[15 15],'court found','BoxColor','green','BoxOpacity',0.4,'TextColor','white');
    end
    %figure ; imshow(tempImage);
    %writeVideo(v2,tempImage);
    %pause(1/vidObj.FrameRate);
end
close(v);
%close(v2);
%use the flag_lst to segment out the video into play and non-play.
toc


end