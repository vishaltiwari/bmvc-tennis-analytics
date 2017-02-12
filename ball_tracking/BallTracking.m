%Write the candidates detection output file.
%video_file = './TRGMCOutputVideo/Test_sample1.avi';
video_file = 'ao_2017_stable_motion_compensation.avi';
candidate_video_out = './TEST_CANDIDATES_2017.avi';
tracklet_video_out = './TEST_TRACKLETS_2017.avi';
ballDetection_video_out = './TEST_BALL_DETECTION_ao_2017_3_R_40_d_th_10.avi';
background_in = '../court_detection/non_court_template_8.png';
path_threshold = 200;

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

N = 10;
R = 40; % 20 -> 30 -> 45 
m_th = 6;
d_th = 5; % 5 to 10 -> 20
k_th = 20;
alpha = 100;

disp('Generating Ball candidates...');
[candidate_per_frame , mov , len, buffer_len] = temporalDifference_mark3(video_file,candidate_video_out,background_in);
disp('Done Generating Ball candidates!!!');

disp('Creating Tracklets by optimizing motion models...');
node_set = CreateTracklet(candidate_per_frame,mov,len,buffer_len,tracklet_video_out,N , R , m_th , d_th);
disp('Done Creating Tracklets!!!');

disp('Creating graph from generated Models...');
graph = createGraphFromNodeList(node_set , k_th);
disp('Done Creating graph!!!');

disp('Getting all Possible Paths from the Graph...');
[allPaths , allWeights] = getAllPairPaths(graph);

disp('Getting the length of all possible paths...');
%[allLengths , node_set_withL] = getLengths(allPaths , node_set);
allLengths = getAllLengthPaths(allPaths , node_set);
disp('Reterived all Possible paths!!!');

disp('Sorting all Paths...');
[SP , SW , SL] = Quicksort_Paths_Weight_Length(allPaths , allWeights, allLengths , path_threshold , alpha);

disp('Get the best paths based on compatibility of paths...');
[bestP2 , bestW , bestL] = getFinalTracks_APSP(SP , SW , SL , node_set);

disp('Writing the detected ball position to a Video file...');
ballPositions = getAllBallPositions(node_set , mov , bestP2 , ballDetection_video_out);

disp('Done!!!');
