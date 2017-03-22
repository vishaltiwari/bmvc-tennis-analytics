run('./vlfeat-0.9.20/toolbox/vl_setup');

%% Test Train Split
v = VideoReader('Federer-Isner-QtrFinals.mp4');
num_frames = round(v.Duration*v.FrameRate);
stride_frames = 5;

%% Reading the label file
%label_file = fopen('Federer-Isner-QtrFinals-hh.txt');
%C = textscan(label_file,'%s\t%s\t%s\t%s\t%d');
%[~, ~, ~, H, MN, S] = datevec( C{2} );
%begin_frames = uint32(round((H*3600 + MN*60 + S)*v.FrameRate));
%[~, ~, ~, H, MN, S] = datevec( C{3} );
%end_frames = uint32(round((H*3600 + MN*60 + S)*v.FrameRate));

label_file = fopen('Federer-Isner-QtrFinals.txt');
C = textscan(label_file,'%s\t%d\t%d\t%d\t%d');
begin_frames = uint32(round((C{2}/1000)*v.FrameRate));
end_frames = uint32(round((C{3}/1000)*v.FrameRate));
play_frames = [ begin_frames, end_frames ];

fclose(label_file);

%% Convert the label file into a vector.
vid_label = zeros(num_frames, 1);
for i = 1:size(play_frames,1)
    vid_label(play_frames(i, 1):play_frames(i, 2)) = 1; % Mark 1 for each play segment
end

% v = VideoReader('Federer-Isner-QtrFinals.mp4');
% i = 0;
% while hasFrame(v)
%     i = i + 1
%     frame = readFrame(v);
%     if i > num_frames
%         break;
%     end
%     if mod(i, stride_frames) ~= 0
%        continue
%     end
%     idx = i/stride_frames;
%     if vid_label(i) == 1
%         imshow(frame);
%     end
% end
% return;

%% Reading the video file, extracting features and saving the features with label
data   = zeros(6820, ceil(num_frames/stride_frames));
label  = zeros(ceil(num_frames/stride_frames), 1);
i = 0;
while hasFrame(v)
    i = i + 1
    frame = readFrame(v);
    if i > num_frames
        break;
    end
    if mod(i, stride_frames) ~= 0
       continue
    end
    idx = i/stride_frames;
    frame = imresize(frame, 0.25);
    fr_hog = vl_hog(single(frame), 16);
    fr_hog = fr_hog(:);
    data(:,idx) = fr_hog; % the feature from the frame
    label(idx) = vid_label(i); % label, coz we are striding 
end

%% Train and Test Split - 50/50
train_frames = ceil(num_frames/(stride_frames*2));
train_data   = data(:,1:train_frames);
train_label  = label(1:train_frames); 
test_data    = data(:,train_frames+1:end);
test_label   = label(train_frames+1:end);

save('segmentdata.mat','train_data','train_label', 'test_data', 'test_label','-v7.3');
