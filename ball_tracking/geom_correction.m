%% TODO: Before running this code, follow:
%   - addpath('/home/vishal/libraries/ransac_homography','-end')
%   - run('../../vlfeat-0.9.20/toolbox/vl_setup')
%% Get the geometric transformation:

%video_file = '/media/vishal/0804E66104E6516C/Racket_games/Lawn_Tennis/semi_sample.mp4'; % Lawn Tennis
% use the de-interlaced video:

clear;

video_file = 'semi_sample_deinterlaced.avi'; % Lawn Tennis

vidObj = VideoReader(video_file);

vidHeight = vidObj.Height;
vidWidth = vidObj.Width;

flag = 0;
prev_frame = [];
prev_frame_feature = [];
prev_frame_descriptor = [];

while hasFrame(vidObj)
    
    % As of now will be using the prev_frame_features as the base onto
    % which to project.
    if flag == 0
        prev_frame = readFrame(vidObj);
        
        prev_frame = single(rgb2gray(prev_frame));
        
        [prev_frame_feature, prev_frame_descriptor] = vl_sift(prev_frame , 'PeakThresh' , 8);
        flag = 1;
        continue;
    end
    
    curr_frame = readFrame(vidObj);
    curr_frame = single(rgb2gray(curr_frame));
    [curr_image_feature, curr_image_descriptor] = vl_sift(curr_frame , 'PeakThresh' , 8);
    
    %% Get the matching features
    [matches, scores] = vl_ubcmatch(curr_image_descriptor, prev_frame_descriptor);
    
    curr_image_feature_matches_index = matches(1,:);
    base_image_feature_matches_index = matches(2,:);
    
    curr_image_features_match = curr_image_feature(1:2,curr_image_feature_matches_index);
    %curr_image_features_match = [curr_image_features_match ; ones(1,size(curr_image_features_match,2))];
    
    base_image_features_match = prev_frame_feature(1:2,base_image_feature_matches_index);
    %base_image_features_match = [base_image_features_match ; ones(1,size(base_image_features_match,2))];
    
    %prev_frame_feature = curr_image_feature;
    %prev_frame_descriptor = curr_image_descriptor;
    
    %% Get the homograpy matrix:
    
    % Use the transformation type:
    tform = fitgeotrans(curr_image_features_match', base_image_features_match' , 'projective');
    
    H = vgg_H_from_x_lin(curr_image_features_match, base_image_features_match);
    H = transpose(H);
    tform = projective2d(H);
    
    %[ H , wim1 ] = homography( curr_frame, prev_frame);
    
    %[H , inlinerPoints] = findHomography(curr_image_coords, base_image_coords);
    
    %imshow(wim1,[]);
    % get the x,y coordinates of the points and plot:
    
    %% Transform the image using the homograph:
    
    rectified_image = imwarp(curr_frame,tform);
    
    %imshow(curr_frame,[]);
    imshow(rectified_image,[]);
    pause(1/vidObj.FrameRate);
end