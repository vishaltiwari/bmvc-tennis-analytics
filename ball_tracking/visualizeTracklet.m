function [frame_image] = visualizeTracklet(frame_image, model,candidates_per_frame,k_start,k_end)
    a = model.compute_a();
    v1 = model.compute_v1(a);
    track = [];
    
    for i=k_start:k_end
        cand = candidates_per_frame{i};
        if ~isempty(cand)
            for j=1:size(cand,2)
                point = cand{1};
                frame_image = insertShape(frame_image,'circle',[point.x point.y 2]);
            end
        end
    end
    
    % extra-polating the values of the model, and then plot them
    k_min = model.p1.k;
    k_max = model.p3.k;
    for i=k_min:k_max
        %frameDist = i - curr_k;
        p_dash = model.compute_k_prime(v1,a,i);
        track = [track ; p_dash];
    end
        %Just plot in the bounds:
%         for i=model.p1.k:model.p3.k
%             %frameDist = i - curr_k;
%             p_dash = model.compute_k_prime(v1,a,i);
%             track = [track ; p_dash];
%         end

    % convert the track into the M matrix.
    M = [];
    for i=1:size(track,1)
        point = track(i,:);
        M = [M point];
    end
    
    frame_image = insertShape(frame_image , 'Line' , M , 'LineWidth' , 5);
    point1 = model.p1;
    point2 = model.p2;
    point3 = model.p3;
    frame_image = insertShape(frame_image, 'circle' , [point1.x point1.y 5] , 'Color' , 'green' , 'Opacity',1);
    frame_image = insertShape(frame_image, 'circle' , [point2.x point2.y 5] , 'Color' , 'green' , 'Opacity',1);
    frame_image = insertShape(frame_image, 'circle' , [point3.x point3.y 5] , 'Color' , 'green' , 'Opacity',1);
    %scatter(track(:,1) , track(:,2) , 100, '.');
    %imshow(frame_image);
%     hold on;
    %plot(track(:,1) , track(:,2),'color',[0 0 0]);
%     scatter(track(:,1) , track(:,2));
%     plot(track(:,1) , track(:,2),'color',[0 0 0]);
    
    % plot the start and end point
    %scatter(track(1,1) , track(1,2) , 10 , 'd');
    %scatter(track(size(track,1),1) , track(size(track,1),2) , 's');
    
end