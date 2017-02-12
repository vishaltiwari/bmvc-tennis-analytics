function [M] = drawTracklet(model,start_k,end_k , flag_track)
    a = model.compute_a();
    v1 = model.compute_v1(a);
    track = [];
    
    % extra-polating the values of the model, and then plot them
    if flag_track == 0
        for i=start_k:end_k
            %frameDist = i - curr_k;
            p_dash = model.compute_k_prime(v1,a,i);
            track = [track ; p_dash];
        end
    else
        %Just plot in the bounds:
        for i=model.p1.k:model.p3.k
            %frameDist = i - curr_k;
            p_dash = model.compute_k_prime(v1,a,i);
            track = [track ; p_dash];
        end
    end

    M = [];
    for i=1:size(track,1)
        M = [M track(i,1) track(i,2)];
    end
    
    %scatter(track(:,1) , track(:,2) , 100, '.');
    %plot(track(:,1) , track(:,2),'color',[0 0 0]);
    
    % plot the start and end point
    %scatter(track(1,1) , track(1,2) , 10 , 'd');
    %scatter(track(size(track,1),1) , track(size(track,1),2) , 's');
    
end