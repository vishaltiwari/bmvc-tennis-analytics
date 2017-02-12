function [X , Y] = extractFeatures(binary_image, curr_image)
    
    CC = bwconncomp(binary_image);
    stats = regionprops(CC,'centroid' , 'Eccentricity' , 'MajorAxisLength' , 'MinorAxisLength', 'Solidity' , 'BoundingBox');
    components = CC.PixelIdxList;
    centroids = cat(1,stats.Centroid);
    majorAxisLength = cat(1,stats.MajorAxisLength);
    minorAxisLength = cat(1,stats.MinorAxisLength);
    bbox = cat(1,stats.BoundingBox);
    solidity = cat(1,stats.Solidity);
    hsv_image = rgb2hsv(curr_image); 
    
    % Extract H,S,V colors:
    Hue = hsv_image(:,:,1);
    Saturation = hsv_image(:,:,2);
    Intensity = hsv_image(:,:,3);

    
    image_bb = insertShape(curr_image,'Rectangle',bbox,'linewidth',1);
    imshow(image_bb);

    X = [];
    Y = [];
    s = CC.ImageSize;
    x = [];
    i = 1;
    while i<=size(components,2)
%    for i=1:size(components,2)
        image_bb = insertShape(curr_image,'Rectangle',bbox,'linewidth',2);
        pixelGroup = components(i);
        allPixelsIndx = pixelGroup{1};
        groupSize = size(allPixelsIndx,1);
        %draw a bounding box around the blob:
        if groupSize < 5
            i = i+1;
            continue
        end
        
        %[rows , cols] = ind2sub(s,allPixelsIndx);
        % color the box red
        curr_box = bbox(i,:);
        image_bb = insertShape(image_bb,'Rectangle',curr_box,'linewidth',2,'Color','red');
        imshow(image_bb);
        
        % Ask the user, if to consider the box to extract features.
        in = input('Do you want to extract features for this box: (1: Yes , 0:No , 2: redo last box , 3: exit funtion, return saved fetures)');
        if in == 0
            i = i + 1;
            continue
        elseif in==2
            X(i-1,:) = [];
            Y(i-1,:) = [];
            i = i-1;
        elseif in==3
            return;
        elseif in~=0 && in~=1 && in~=2 && in~=3 
            continue;
        end
        
        in = input('Is it a ball or a non-ball: (1: Yes , 0: No)');
        out = 0;
        if in == 1
            out = 1;
        else
            out = 0;
        end
        
        center = centroids(i,:);
        majorAxis = majorAxisLength(i);
        minorAxis = minorAxisLength(i);
        solid = solidity(i);
        
        mean_H = mean(Hue(allPixelsIndx));
        mean_S = mean(Saturation(allPixelsIndx));
        mean_I = mean(Intensity(allPixelsIndx));
        
        x = [mean_H , mean_S , mean_I , center , majorAxis , minorAxis , solid];
        X = [X ; x];
        Y = [Y ; out];
        % change it back to yellow:
        image_bb = insertShape(image_bb,'Rectangle',curr_box,'linewidth',1,'Color','yellow');
        i = i + 1;
    end
    
end