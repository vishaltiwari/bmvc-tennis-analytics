% get detection accuracy:

seq_list = [2,3,5,7,9,19,20];
total_det = 0;
total_frames = 0;
for blah=1:size(seq_list,2)
seq_nb = seq_list(blah);
%seq_nb = 9;
filename_test = ['test_data/',num2str(seq_nb),'_ballpoints.mat'];
load(filename_test);
filename_ballPositions = ['result/' , num2str(seq_nb) , '_ballPositions.mat'];
load(filename_ballPositions);
buffer_len=6;
non_det = 0;
det = 1;
ptr = 1;
vid = VideoReader('result/2_TEST_BALL_DETECTION.avi');
frame = readFrame(vid);
% figure ; imshow(frame);
% hold on;
for i=1:size(ball_points,1)
    actual_point = ball_points(i,:);
    frame_nb = actual_point(3) - buffer_len;
    [pred_point]= getBallPositionAtFrame(ballPositions,frame_nb,buffer_len);
    scatter(pred_point.x , pred_point.y , '+');
    scatter(actual_point(1) , actual_point(2) , 'o');
    %pred_point = ballPositions(frame_nb-17+1);
    dis = sqrt((pred_point.x - actual_point(1))^2 + (pred_point.y - actual_point(2)));
    if dis>10
        non_det = non_det + 1;
    else
        det = det + 1;
    end
end
total_det = total_det + det;
total_frames = total_frames + size(ball_points,1);
disp(['Seq no: ' , num2str(seq_nb) , '  total_frames:' , num2str(size(ball_points,1)) , ' accuracy: ', num2str(det / size(ball_points,1))]);
end
disp(['total frames: ' , num2str(total_frames)]);
disp(['Overall Accuracy:' , num2str(total_det / total_frames)]);