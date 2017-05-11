function [candidate_per_frame , mov , len, buffer_len] = temporalDifference_mark3(video_file,video_out,background_in,write_flag,upper_player_bbox , result_frame_upper , lower_player_bbox , result_frame_lower)
%video_file = 'ov_sample1.avi'; % sample1.avi
%video_file = 'ov_sample2.avi'; % sample2.avi
%video_file = './TRGMCOutputVideo/ov_aus_open_sample1.avi'; % AUS OPEN
%video_file = './TRGMCOutputVideo/australian_open_sample2.avi'; % wimbelden

% Write to an output files:
%video_out = './TRGMCOutputVideo/australian_open_sample2_CANDIDATES.avi';

%% load the player_bbox from ../player_tracking_color_seg/player_bbox.mat
%load('../player_tracking_color_seg/upper_lower_player_bbox.mat');
% variables are upper_player_bbox and lower_player_bbox.

background = imread(background_in);
background = imbinarize(rgb2gray(background));

vidObj = VideoReader(video_file);
vidOutObj = VideoWriter(video_out);
open(vidOutObj);

vidHeight = vidObj.Height;
vidWidth = vidObj.Width;

buffer_len = 6; % 6 frames on both the side current frame

% Get the first 2 * buffer_len video frame array (sliding window over time)
%import java.util.LinkedList;
%time_window_buffer = LinkedList();

% Reading all the frames and storing them all

mov = struct('cdata',zeros(vidHeight,vidWidth,3,'uint8'),...
    'colormap',[]);

len=1;
% pad the frames with ones on both ends:
for a=1:buffer_len
    ones_frame = ones(vidHeight , vidWidth,3);
    mov(len).cdata = ones_frame;
    len = len + 1;
end
frame_count = 0;
while hasFrame(vidObj)
    mov(len).cdata = readFrame(vidObj);
    len = len+1;
    frame_count = frame_count + 1;
end

for a=1:buffer_len
    ones_frame = ones(vidHeight , vidWidth,3);
    mov(len).cdata = ones_frame;
    len = len + 1;
end

start_ptr = 1;
curr_ptr = buffer_len + 1;
end_ptr = 2*buffer_len + 1;
total_frame = len;

threshold = 10; % NOTE: A hyperparamter

final_features = [];
final_label = [];
while (curr_ptr-buffer_len) < frame_count %&& (curr_ptr-buffer_len) < size(upper_player_bbox,1)
    center_image = imgaussfilt(rgb2gray(mov(curr_ptr).cdata),1,'FilterSize',5);
    binary_image = ones(size(center_image));
    for j=-6:2:6
        if j == -2 || j ==0 || j == 2
            continue
        end
        subtracting_image = imgaussfilt(rgb2gray(mov(curr_ptr + j).cdata),1,'FilterSize',5);
        diff_image = (center_image - uint8(subtracting_image)) > threshold;
        binary_image = diff_image & binary_image;
    end
    %imshow(binary_image);
    % apply dilation on the image. Erosion then dilation.
    %binary_image = binary_image & background; % Remove the part which is noncenterd.
    binary_image = uint8(binary_image .* 255);
    %binary_image = binary_image & background;
%     IM2 = imerode(binary_image,ones(3,3));
%     IM2 = imdilate(IM2 , ones(3,3));
    %binary_image = uint8(IM2 .* 255);
    index_upper = find(result_frame_upper == (curr_ptr - buffer_len));
    index_lower = find(result_frame_lower == (curr_ptr - buffer_len));
    
    low_box = lower_player_bbox(index_lower,:);
    up_box = upper_player_bbox(index_upper,:);
    candidate_per_frame{curr_ptr - buffer_len} = getBallCandidates_mark3(binary_image,mov(curr_ptr).cdata,low_box , up_box ,(curr_ptr - buffer_len));
    
    %Plot the candidates as small circles on the video:
    if write_flag == 1
        color_image = mov(curr_ptr).cdata;
        cand_lst = candidate_per_frame{curr_ptr - buffer_len};
        for c_indx=1:size(cand_lst,2)
            point = cand_lst{c_indx};
            if isempty(point)
                continue;
            end
            color_image = insertShape(color_image,'circle',[floor(point.x) , floor(point.y) , 5] );
        end
        
        writeVideo(vidOutObj,color_image);
    end
    
    start_ptr = start_ptr + 1;
    curr_ptr = curr_ptr + 1;
    end_ptr = end_ptr + 1;
end
close(vidOutObj);

end