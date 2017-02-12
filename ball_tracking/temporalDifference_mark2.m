function [final_features , final_label] = temporalDifference_mark2(video_file,video_out , background_in)
%video_file = 'ov_sample1.avi'; % sample1.avi
%video_file = 'ov_sample2.avi'; % sample2.avi
%video_file = './TRGMCOutputVideo/ov_aus_open_sample1.avi'; % AUS OPEN
%video_file = './TRGMCOutputVideo/australian_open_sample2.avi'; % wimbelden

% Write to an output files:
%video_out = './TRGMCOutputVideo/australian_open_sample2_CANDIDATES.avi';

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

while hasFrame(vidObj)
    mov(len).cdata = readFrame(vidObj);
    len = len+1;
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

diff_count = 0;
final_features = [];
final_label = [];
while end_ptr < total_frame
    
    if curr_ptr <= 227
        diff_count = diff_count + 1;
        start_ptr = start_ptr + 1;
        curr_ptr = curr_ptr + 1;
        end_ptr = end_ptr + 1;
        continue
    end
    
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
    binary_image = binary_image & background; % Remove the part which is noncenterd.
    binary_image = uint8(binary_image .* 255);
    IM2 = imerode(binary_image,ones(5,5));
    IM2 = imdilate(IM2 , ones(5,5));
    binary_image = uint8(IM2 .* 255);
    [frameFeaturesX , featureLabel]= extractFeatures(binary_image , mov(curr_ptr).cdata);
    final_features = [final_features ; frameFeaturesX];
    final_label = [final_label ; featureLabel];
    writeVideo(vidOutObj,binary_image);
    
    if mod(diff_count,20) == 0
        save('../candidate_model/features_so_far.mat','final_features','final_label','curr_ptr');
    end
    
%    [height , width] = size(IM2);
%     binary_image_rgb = zeros(height,width,3);
%     binary_image_rgb(:,:,1) = uint8(IM2);
%     binary_image_rgb(:,:,2) = uint8(IM2);
%     binary_image_rgb(:,:,3) = uint8(IM2);
%    binary_image_rgb = [IM2; IM2; IM2];
%     final_image = mov(curr_ptr + j).cdata + uint8(binary_image_rgb);
%     % Make an or of prev, curr, and next binray frame.
%     if diff_count == 0
%         prev_diff = IM2;  
%     elseif diff_count == 1
%         curr_diff = IM2;
%         writeVideo(vidOutObj,prev_diff);
%     elseif diff_count >=2
%         next_diff = IM2;
%         combined_image = prev_diff | curr_diff | next_diff;
%         combined_image = uint8(combined_image .* 255);
%         combined_image = imdilate(combined_image , ones(9,9));
%        % morphed_image = 0.8 * mov(curr_ptr-1).cdata + combined_image;
%         s = regionprops(combined_image,'centroid');
%         centroids = cat(1, s.Centroid);
%         %combined_image = insertShape(combined_image,'FilledCircle',[centroids(:,1),centroids(:,2),5]);
%         [idx,C] = kmeans(centroids,3);
%         writeVideo(vidOutObj,combined_image);
%         prev_diff = curr_diff;
%         curr_diff = next_diff;
%     end
%    
%     %writeVideo(vidOutObj,IM2);
%     diff_count = diff_count + 1;
    diff_count = diff_count + 1;
    start_ptr = start_ptr + 1;
    curr_ptr = curr_ptr + 1;
    end_ptr = end_ptr + 1;
end
close(vidOutObj);

end