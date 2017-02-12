function [] = plotTrackletOnImage(coords,curr_image, model,start_k,end_k , flag_track)
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

    %plot all the candidates of the window:
    for i=1:size(coords,1)
        curr_image = insertShape(curr_image, 'circle' ,[coords(i,1) , coords(i,2) , 1] );
    end

    imshow(curr_image);
    %hold on;
    %plot(track(:,1) , track(:,2) , 'color' , [1 1 1]);
    
end