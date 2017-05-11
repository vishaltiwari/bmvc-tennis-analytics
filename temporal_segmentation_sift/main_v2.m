% main2
%segment out the frames:
run('../../vlfeat-0.9.20/toolbox/vl_setup');
%% Input, output, template parameters
filename = 'Federer-vs-Nadal-AO2017';
%filename = 'test_sample4';
video_in = ['sample_data/' , filename , '.mp4'];
%video_write_prefix = ['/media/vishal/0804E66104E6516C/Racket_games/temporal_segments/' , filename];
video_write_prefix = ['/Neutron9/anurag/AO_2017_segments/' , filename];
template_in = 'court_template/template_new.png';
mkdir(video_write_prefix);

vidObj = VideoReader(video_in);
start_time = 360; %in secs
remaing_time = 1440;
end_time = vidObj.Duration - remaing_time;
total_frames = ceil(vidObj.Duration*vidObj.FrameRate);

%Get the features of the template provided.
[templateFeatures , templateDescriptors , court_pos] = getTemplateFeaturesDescriptors(template_in);
base_image = imread(template_in);
base_image = imresize(base_image , 0.10);

tic
prev_frame_court = 0;
i = 0;
seq_arr = {};
not_found_court = 0;
thresh = 10;
seq_to_write = 0;
first_flag = 1;
found_court_len = 0;
min_court_len_thresh = 50;
curr_seq_frames = {};
end_frame_nb = 0;

count = 1;
video_write_avi = [video_write_prefix , '/' , num2str(count) , '.avi'];
video_write_mp4 = [video_write_prefix , '/' , num2str(count) , '.mp4'];
vidOut = VideoWriter(video_write_avi , 'Motion JPEG AVI');
%vidOut = VideoWriter(video_write , 'MPEG-4');
open(vidOut);
c = 0;
while hasFrame(vidObj)
    i = i + 1;
    curr_image_orig = readFrame(vidObj);
    curr_image = imresize(curr_image_orig , 0.10);
    if i > ceil(end_time*vidObj.FrameRate)
        break;
    end
    if i < ceil(start_time*vidObj.FrameRate)
        continue;
    end
    %Get the features of the current frame.
    [frameFeatures , frameDescriptors] = getFrameFeatureDescriptor(curr_image);
    
    %Match template desciptors to the current frame descriptoes.
    [court_flag , tform] = matchDescriptors(templateFeatures , templateDescriptors , frameFeatures , frameDescriptors,curr_image , template_in);

    if court_flag == 0
        if seq_to_write == 1 && not_found_court >= thresh
            if found_court_len > min_court_len_thresh
                for j=1:found_court_len
                    writeVideo(vidOut,curr_seq_frames{j});
                end
                count = count + 1;
                %release(videoFWriter);
                close(vidOut);
                exe_str = ['ffmpeg -i ' ,video_write_avi, ' -c:v libx264 -crf 19 -preset slow -c:a libfdk_aac -b:a 192k -ac 2 ' , video_write_mp4]';
                system(exe_str');
                rem_avi = ['rm -rf ' , video_write_avi];
                system(rem_avi);
                video_write_avi = [video_write_prefix , '/' , num2str(count) , '.avi'];
                video_write_mp4 = [video_write_prefix , '/' , num2str(count) , '.mp4'];
                %videoFWriter = vision.VideoFileWriter(video_write);
                vidOut = VideoWriter(video_write_avi,'Motion JPEG AVI');
                %vidOut = VideoWriter(video_write_avi , 'MPEG-4');
                open(vidOut);
                seq_to_write = 0;
                found_court_len = 0;
            end
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
            found_court_len = 0;
        else
            % add the missing frames.
        end
        prev_frame_court = 1;
        end_frame_nb = i;
        not_found_court = 0;
        found_court_len = found_court_len + 1;
        %step(videoFWriter , curr_image);
        %im_to_write = imresize(curr_image_orig,0.6);
        %curr_seq_frames{found_court_len} = im_to_write;
        curr_seq_frames{found_court_len} = curr_image_orig;
        %writeVideo(vidOut,curr_image_orig);
    end
    
end
toc
close(vidObj);