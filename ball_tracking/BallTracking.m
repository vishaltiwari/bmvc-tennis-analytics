tic
%Write the candidates detection output file.
%video_file = './TRGMCOutputVideo/Test_sample1.avi';
addpath('/Neutron9/anurag/AO_2017_segments/Federer-vs-Nadal-AO2017-stab');
addpath('../player_seg_color/player_detection_results/');
addpath('../player_seg_color/');

seq_nb = 126;
video_file = [num2str(seq_nb),'.avi.stab.avi'];
%video_file = 'ao_2017_stable_motion_compensation.avi';
%video_file = 'test.avi';
candidate_video_out = ['result/',num2str(seq_nb),'_TEST_CANDIDATES.avi'];
tracklet_video_out = ['result/',num2str(seq_nb),'_TEST_TRACKLETS.avi'];

ballDetection_video_out = ['result/',num2str(seq_nb),'_TEST_BALL_DETECTION.avi'];
background_in = '../court_detection/non_court_template_8.png';
path_threshold = 200;

upper_player_bbox_file = ['upper_player_bbox_',num2str(seq_nb),'.txt'];
lower_player_bbox_file = ['lower_player_bbox_',num2str(seq_nb),'.txt'];

result_data_dir = '../player_seg_color/player_detection_results';
player_type = 'upper_player_bbox';
[result_bbox_upper , result_frame_upper] = parsedatafile(seq_nb , result_data_dir , player_type);

player_type = 'lower_player_bbox';
[result_bbox_lower , result_frame_lower] = parsedatafile(seq_nb , result_data_dir , player_type);

result_bbox_lower(:,1) = result_bbox_lower(:,1) - (result_bbox_lower(:,3)*0.5);
result_bbox_lower(:,2) = result_bbox_lower(:,2) - (result_bbox_lower(:,4)*0.5);

result_bbox_upper(:,1) = result_bbox_upper(:,1) - (result_bbox_upper(:,3)*0.5);
result_bbox_upper(:,2) = result_bbox_upper(:,2) - (result_bbox_upper(:,4)*0.5);

%Need to set other hyderparameters here,
% N (size of the sliding window when creating tracklets) (10)
% R (Radius where to look for seed triplets),
% m_th: (motion model optimization), least number of elements in the support
% set.
% buffer_len: (size of the sliding window),
% d_th: (distance to consider for including support point),
% k_th: (graph creation: Max frame over which to join nodes),
% Lth (filtering out the best Paths [not used in code, but used in paper])
% alpha(when comparing two paths).

N = 15;
R = 50; % 20 -> 30 -> 45 -> 60 ->80 ->120
m_th = 4;
d_th = 5; % 5 to 10 -> 20
k_th = 15; % 30 > 40
alpha = 1000000;

disp('Generating Ball candidates...');
[candidate_per_frame , mov , len, buffer_len] = temporalDifference_mark3(video_file,candidate_video_out,background_in,1,result_bbox_upper , result_frame_upper , result_bbox_lower , result_frame_lower);
disp('Done Generating Ball candidates!!!');

disp('Creating Tracklets by optimizing motion models...');
node_set = CreateTracklet(candidate_per_frame,mov,len,buffer_len,tracklet_video_out,N , R , m_th , d_th,0);
disp('Done Creating Tracklets!!!');

disp('Creating graph from generated Models...');
graph = createGraphFromNodeList(node_set , k_th);
disp('Done Creating graph!!!');

disp('Getting all Possible Paths from the Graph...');
[allPaths , allWeights] = getAllPairPaths(graph);

disp('Getting the length of all possible paths...');
[allLengths , node_set_withL] = getLengths(allPaths , node_set);
%allLengths = getAllLengthPaths(allPaths , node_set);
disp('Reterived all Possible paths!!!');

disp('Sorting all Paths...');
[SP , SW , SL] = Quicksort_Paths_Weight_Length(allPaths , allWeights, allLengths , path_threshold , alpha);

disp('Get the best paths based on compatibility of paths...');

%[bestP2 , bestW , bestL] = getFinalTracks_APSP(SP , SW , SL , node_set_withL);
bestP2 = SP(1);
bestW = SL(1);
bestL = SW(1);

% draw the bestP2
disp('Getting the ball positions and showing it on video...');
ballPositions = writeBestPath(node_set , mov , len , buffer_len , bestP2 , ballDetection_video_out)

save(['result/',num2str(seq_nb) , '_ballPositions.mat'] , 'ballPositions');
%disp('Writing the detected ball position to a Video file...');
%ballPositions = getAllBallPositions(node_set , mov , bestP2 , ballDetection_video_out);


disp('Done!!!');
toc