run('./vlfeat-0.9.20/toolbox/vl_setup');

video_files = [ {'Baghdatis-Murray-R16'}, {'Federer-Delpotro-SemiFinals'}, ...
                {'Hewitt-Djokovic-R16'}, {'Murray-Djokovic-SemiFinals'}, ...
                {'Murray-Federer-Finals'} ];
            
%for fname=video_files
%fname = fname{:}

%video_path = ['/home/anurag/TennisLondon/', fname, '.mp4'];
%model_path = 'SVM-C0.05KChi2_O2_model.mat';
%rally_path = ['./rally/', fname, '/'];

% load(model_path);
% i = 0;
% v = VideoReader(video_path);
% num_frames = round(v.Duration*v.FrameRate);
% vid_label = zeros(num_frames, 1);
% hom.kernel = 'KChi2';
% hom.order = 2;
% while hasFrame(v)
%     i = i + 1;
%     if mod(i, 1000) == 0
%         disp(['Iter :', num2str(i)]);
%     end
%     frame = readFrame(v);
%     frame = imresize(frame, 0.25);
%     fr_hog = vl_hog(single(frame), 16);
%     fr_hog = fr_hog(:);
%     val = vl_svmdataset(fr_hog, 'homkermap', hom);
%     [~,~,~, scores] = vl_svmtrain(val, [-1], 0, 'model', w, 'bias', b, 'solver', 'none') ;
%     scores(scores > 0) = 1;
%     scores(scores < 0) = -1;
%     vid_label(i) = scores(1);
% end
% 
% save([ fname, '.mat'],'vid_label', '-v7.3');

%end

fname = video_files{5};
video_path = ['/home/anurag/TennisLondon/', fname, '.mp4'];
model_path = 'SVM-C0.05KChi2_O2_model.mat';
rally_path = ['./rally/', fname, '/'];
load([ fname, '.mat']);

v = VideoReader(video_path);

i = 0;
seg_frames = 0;
num_vid = 0;
flag = 0;
frame_set = {};
while hasFrame(v)
    i = i + 1;
    frame = readFrame(v);
    if vid_label(i) == 1
        seg_frames = seg_frames + 1;
        if seg_frames == 1
            flag = 1;
        end
        frame_set{seg_frames} = frame;
    end
    if vid_label(i) == -1 && flag == 1
        disp([ 'Rally ', num2str(num_vid) ]);
        flag = 0;
        seg_frames = 0;
        num_vid = num_vid + 1;
        vid = VideoWriter([ rally_path, 'rally-', num2str(num_vid), '.mp4']);
        vid.FrameRate = 25;
        open(vid);
        for frame=frame_set
            writeVideo(vid, frame{:});
        end
        close(vid);
        clearvars frame_set
        frame_set = {};
    end
end