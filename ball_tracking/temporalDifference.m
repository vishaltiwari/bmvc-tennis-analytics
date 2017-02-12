function [candidate_per_frame , mov , len, buffer_len] = temporalDifference(video_file,video_out,background_in)

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

time_window_buffer = [];
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

i = buffer_len+1;
candidate_per_frame = cell(1,len);
count = 1;

%% TODO:
% TODO: Add 7 frames, where every pixels is a '1',
% on both side of the time scales,
% so as to retain the output video of the same size.

for k=1:len-1
    frame_ahead = imgaussfilt(rgb2gray(mov(k).cdata),1,'FilterSize',5);
    if k <= 2*buffer_len+1
        time_window_buffer{k} = frame_ahead;
        continue;
    end
    frame_curr = imgaussfilt(rgb2gray(mov(i).cdata),1,'FilterSize',5);
    
    a = 1;
    diff_collection = [];
    binary_image = ones(size(frame_curr));
    threshold = 10;
    for j=-6:2:6
        if j == 0 || j == 2 || j == -2
            continue;
        end
        diff = (frame_curr - uint8(time_window_buffer{j+buffer_len+1})) > threshold;
        binary_image = diff & binary_image;
        a = a + 1;
    end
    
    for j=1:2* buffer_len
        time_window_buffer{j} = time_window_buffer{j+1};
    end
    time_window_buffer{2*buffer_len+1} = frame_ahead;
    
    % move one the time filter one step
    % time_window_buffer.remove()
    
    %imshow(binary_image);
    
    binary_image_boundry = bwmorph(binary_image,'remove');
    %binary_image_boundry = edge(binary_image,'sobel');
    
    %imshow(binary_image_boundry);
    %% Get the ball candidate from temporal difference image.
    %Get the x,y gradient:
    %[Gmag,Gdir] = imgradient(frame_curr);
    Gdir = [];
    %imshow(binary_image_boundry);

    %candidate_per_frame{count} = ballCandidate(binary_image_boundry,vidHeight,count);
    candidate_per_frame{count} = getBallCandidates(binary_image,mov(i).cdata , count);
    
    %Plot the candidates as small circles on the video:
    color_image = mov(i).cdata;
    cand_lst = candidate_per_frame{count};
    for c_indx=1:size(cand_lst,2)
        point = cand_lst{c_indx};
        if isempty(point)
            continue;
        end
        color_image = insertShape(color_image,'circle',[floor(point.x) , floor(point.y) , 5] );
    end
%    imshow(color_image);
    writeVideo(vidOutObj,color_image);
    pause(1/vidObj.FrameRate);
    
    i = i+1;
    count = count + 1;
end

close(vidOutObj);

end