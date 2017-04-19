%function [] = videoGlobalMotionCompensation(video_in , video_write , template_in)
%Need to load/setup vl_feat before running this.
run('../../vlfeat-0.9.20/toolbox/vl_setup');
%% Input, output, template parameters
filename = 'test_sample3';
video_in = ['test_data/' , filename , '.mp4'];
video_write_prefix = ['temporal_segments/' , filename];
mkdir(video_write_prefix);
template_in = 'court_template/template_new.png';

tic
vidObj = VideoReader(video_in);

%Get the features of the template provided.
[templateFeatures , templateDescriptors , court_pos] = getTemplateFeaturesDescriptors(template_in);
base_image = imread(template_in);
base_image = imresize(base_image , 0.10);

%hdinterlacer = vision.Deinterlacer;
count = 1;
video_write = [video_write_prefix , '/' , num2str(count) , '.avi'];
vidOut = VideoWriter(video_write , 'Motion JPEG AVI');
open(vidOut);
%videoFWriter = vision.VideoFileWriter(video_write);
prev_frame_court = 0;
seq_frame = []; %start_frame_nb end_frame_nb
i = 0; % this is the frame number.
not_found_court = 0;
thresh = 5;
seq_to_write = 0;
first_flag = 1;
found_court_len = 0;
min_court_len_thresh = 10;
curr_seq_frames = {};
while hasFrame(vidObj)
    curr_image_orig = readFrame(vidObj);
    curr_image = imresize(curr_image_orig , 0.10);

    i = i + 1;
    if i < 250
        continue;
    end

    %Get the features of the current frame.
    [frameFeatures , frameDescriptors] = getFrameFeatureDescriptor(curr_image);
    
    %Match template desciptors to the current frame descriptoes.
    [court_flag , tform] = matchDescriptors(templateFeatures , templateDescriptors , frameFeatures , frameDescriptors,curr_image , template_in);
    if court_flag == 0
        if seq_to_write == 1 && not_found_court >= thresh
            if found_court_len < min_court_len_thresh
                
            end
            count = count + 1;
            %release(videoFWriter);
            close(vidOut);
            video_write = [video_write_prefix , '/' , num2str(count) , '.avi'];
            %videoFWriter = vision.VideoFileWriter(video_write);
            vidOut = VideoWriter(video_write,'Motion JPEG AVI');
            open(vidOut);
            seq_to_write = 0;
        end
        not_found_court = not_found_court + 1;
        prev_frame_court = 0;
        continue
    else
        if prev_frame_court == 0 && not_found_court >= thresh || first_flag == 1
            frame_count = 1;
            start_frame_nb = i;
            seq_to_write = 1;
            first_flag = 0;
        end
        prev_frame_court = 1;
        end_frame_nb = i;
        not_found_court = 0;
        found_court_len = found_court_len + 1;
        %step(videoFWriter , curr_image);
        writeVideo(vidOut,curr_image_orig);
        %invtform = invert(tform);
        % deinterlace the frame
        %curr_image = step(hdinterlacer , curr_image);
        %image_new = imwarp(curr_image , invtform , 'OutputView',imref2d(size(base_image)));
        %writeVideo(v,image_new);
    end
end
close(vidObj);
if ~isempty(vidOut)
    close(vidOut);
end

toc


%end