function[candidates] = getBallCandidates_mark3(binary_diff,frame_image,lower_player_bbox, upper_player_bbox,frame_nb)
%     blah = insertShape(double(binary_diff),'rectangle',upper_player_bbox);
%     blah = insertShape(blah,'rectangle',lower_player_bbox);
%     imshow(blah);
    candidates = {};
    [height , width ] = size(binary_diff);
    mask = ones(height , width);
    % add box for lower_player
    % width and height on both the boxes much be non-negative.
    
    if ~isempty(upper_player_bbox) && upper_player_bbox(3) > 0 && upper_player_bbox(4) > 0 
        mask(upper_player_bbox(2):upper_player_bbox(2)+upper_player_bbox(4) , upper_player_bbox(1):upper_player_bbox(1)+upper_player_bbox(3)) = 0;
    end
    if ~isempty(lower_player_bbox) && lower_player_bbox(3) > 0 && lower_player_bbox(4) > 0
        mask(lower_player_bbox(2):lower_player_bbox(2)+lower_player_bbox(4) , lower_player_bbox(1):lower_player_bbox(1)+lower_player_bbox(3)) = 0;
    end
    if sum(size(mask) == size(binary_diff)) ~= 2
        mask = ones(height , width);
    end
    new_binary = binary_diff & mask;
    
    new_binary(1:floor(height*0.10) , :) = 0;
    new_binary(floor(height - height*0.10):height , :) = 0;
    %new_binary = imclose(new_binary,ones(3,3));
    %new_binary = medfilt2(new_binary);
    %imshow(new_binary);
    
    % Get all the connected components:
    
    stats = regionprops(new_binary,'centroid' , 'PixelIdxList');
    
    nblobs = size(stats);
    max_pixel_count = 100;
    min_pixel_count = 5;

    count = 1;
    
    for i=1:nblobs
        pixelGroup = stats(i).PixelIdxList;
        groupSize = size(pixelGroup,1);
        %Filter on size.
        if groupSize >= max_pixel_count || groupSize <=min_pixel_count
            continue;
        end
        
        center = stats(i).Centroid;
        p = Point(floor(center(1)),floor(center(2)),frame_nb);
        
        candidates{count} = p;
        count = count + 1;
    end    
end