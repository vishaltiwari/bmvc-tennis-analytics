function [node_set] = CreateTracklet (candidate_per_frame,mov,len,buffer_len,video_out , N , R, m_th , d_th , writeFlag)
%plot the candidate over time.

%video_out = './TRGMCOutputVideo/ov_wim_sample2.avi_TrackLets.avi';

vidOutObj = VideoWriter(video_out);
open(vidOutObj);

% Write the tracklet_formation
vidOutObjTracklet = VideoWriter('tracklet_formation.avi');
open(vidOutObjTracklet);


if isempty(candidate_per_frame)
    return;
end
%Changing this from 7 to 10.
%N = 10;
% Iterate over all the frames, slide over a window of size 2N+1;
ahead_ptr = 2*N+1;
prev_ptr = 1;
% start with current frame
%tracklets_arr = [];
%count = 1;
model_set = [];
total_models = 0;
node_set = [];

img_ptr_mov = buffer_len + 1; % The position in the mov structure.
tracklet_count = 0;
total_frames = size(candidate_per_frame,2);
for k_curr=N+1:total_frames-1
    
    if ahead_ptr > total_frames
        break;
    end
    
    coords = []; %coords(0)::x , coords(1)::y
    
    for j=prev_ptr:ahead_ptr
        candidates = candidate_per_frame{j};
        for c=1:size(candidates,2)
            if isempty(candidates{c})
                continue;
            end
            coords = [coords ; candidates{c}.x candidates{c}.y candidates{c}.k];
        end
    end
    %handle = scatter(coords(:,1) , coords(:,2) , 'b');
    % From this set of data, create tracklets for frame k_curr:
    curr_frame_candidates = candidate_per_frame{k_curr};
    one_prev_frame_candidates = candidate_per_frame{k_curr-1};
    one_next_frame_candidates = candidate_per_frame{k_curr+1};
    
    %Get all coordinates of one prev frame.
    one_prev_coords = [];
    for c=1:size(one_prev_frame_candidates,2)
        if isempty(one_prev_frame_candidates{c})
            continue;
        end
        one_prev_coords = [one_prev_coords ; one_prev_frame_candidates{c}.x  one_prev_frame_candidates{c}.y];
    end
    
    % Get all coordinates of one frame next.
    one_next_coords = [];
    for c=1:size(one_next_frame_candidates,2)
        if isempty(one_next_frame_candidates{c})
            continue;
        end
        one_next_coords = [one_next_coords ; one_next_frame_candidates{c}.x  one_next_frame_candidates{c}.y];
    end
    
    %make the one prev and one next frames kdTrees
    if isempty(one_next_coords) || isempty(one_prev_coords) 
        continue;
    end
    prevKdTree = KDTreeSearcher(one_prev_coords);
    nextKdTree = KDTreeSearcher(one_next_coords);
    
    curr_coords = [];
    for c=1:size(curr_frame_candidates,2)
        if isempty(curr_frame_candidates{c})
            continue;
        end
        curr_coords = [curr_coords ; curr_frame_candidates{c}.x curr_frame_candidates{c}.y];
    end
    if isempty(curr_coords)
        prev_ptr = prev_ptr + 1;
        ahead_ptr = ahead_ptr + 1;
        disp('is empty');
        continue;
    end
    % Get the coordinates in prev,
    [prev_cand_indx , prev_cand_D] = knnsearch(prevKdTree , curr_coords , 'K' , 1);
    [next_cand_indx , next_cand_D] = knnsearch(nextKdTree , curr_coords , 'K' , 1);
        
    %% Filter by R, and get seed pixels
    
    %Get the seeds, and the initial models
    models_lst = [];
    for i=1:size(curr_coords)
        % Create a model of the three seed points:
        if prev_cand_D(i) <= R && next_cand_D(i) <=R
            p1 = Point(one_prev_coords(prev_cand_indx(i),1) , one_prev_coords(prev_cand_indx(i),2) , k_curr-1);
            p2 = Point(curr_coords(i,1) , curr_coords(i,2) , k_curr);
            p3 = Point(one_next_coords(next_cand_indx(i),1) , one_next_coords(next_cand_indx(i),2) , k_curr+1);
            %note in model class there is no check for p1 < p2 < p3 (in frame domain)
            model = Model(p1,p2,p3);
            model.support_set = [p1,p2,p3];
            models_lst = [models_lst model];
            tracklet_count = tracklet_count + 1;
        end
    end
    
    %% For each model in the model_lst, optimize the model
    optimized_models_lst = optimizeMotionModel(models_lst,candidate_per_frame,prev_ptr , ahead_ptr , m_th , d_th , mov(buffer_len+1).cdata , vidOutObjTracklet , buffer_len , 0);
    total_models = total_models + size(optimized_models_lst,2);
    model_set = [model_set optimized_models_lst];
    %convert the optimized models into nodes:
    for i=1:size(optimized_models_lst,2)
        node = Node(optimized_models_lst(i) , k_curr);
        node_set = [node_set node];
    end
    
    %tracklets_arr{count} = optimized_models_lst;
    %% After the optimization:

    disp(tracklet_count); 
    %display
    if writeFlag
        %plot all the candidates of the window:
        curr_image = mov(k_curr).cdata;
        for i=1:size(coords,1)
            curr_image = insertShape(curr_image, 'circle' ,[coords(i,1) , coords(i,2) , 1] );
        end

        %imshow(curr_image);
        allMs = {};
        count = 1;
        for i=1:size(optimized_models_lst,2)
            %hold on ;
            %visualizeTracklet(optimized_models_lst(i),prev_ptr,ahead_ptr,1);
            M = drawTracklet(optimized_models_lst(i),prev_ptr,ahead_ptr,1);
            allMs{count} = M;
            count = count + 1;
            % show the tracklets on the image.
            %plotTrackletOnImage(coords,mov(img_ptr_mov).cdata , models_lst(i) , prev_ptr , ahead_ptr , 0);
            
        end
        if ~isempty(allMs)
            curr_image = insertShape(curr_image,'Line',allMs,'LineWidth',2,'Color',{'red'});
        end
        writeVideo(vidOutObj,curr_image);
        
    end
    img_ptr_mov = img_ptr_mov + 1;
    
    prev_ptr = prev_ptr + 1;
    ahead_ptr = ahead_ptr + 1;
end

close(vidOutObj);
close(vidOutObjTracklet);
%create the directed graph of tracklets.
%graph = createGraphFromList(model_set);
%use model set:

end