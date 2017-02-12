function [candidates] = getBallCandidates(binary_image , curr_img , k)

    candidates = {};
    count = 1;
    addpath('/home/vishal/libraries/fit_ellipse');
    
    hsv_image = rgb2hsv(curr_img);
    
    % Extract H,S,V colors:
    Hue = hsv_image(:,:,1);
    Saturation = hsv_image(:,:,2);
    Intensity = hsv_image(:,:,3);


    %se = ones(7,7);
    %t = imdilate(binary_image , se);
    %binary_image = imerode(t,se);
    %imshow(binary_image);
    CC = bwconncomp(binary_image);
    components = CC.PixelIdxList;
    
    max_pixel_count = 100;
    
    s = CC.ImageSize;
    for i=1:size(components,2)
        pixelGroup = components(i);
        allPixelsIndx = pixelGroup{1};
        [rows , cols] = ind2sub(s,allPixelsIndx);
        
        groupSize = size(allPixelsIndx,1);
        
        %Filter on size.
        if groupSize >= max_pixel_count
            continue;
        end
        
%         mean_H = mean(Hue(allPixelsIndx));
%         mean_S = mean(Saturation(allPixelsIndx));
%         mean_I = mean(Intensity(allPixelsIndx));
% 
%         % Filter on the values, collected from the ball samples.
%         if mean_H < 0.4 || mean_H > 0.6
%             continue
%         end
%         if mean_S < 0.28 || mean_S > 0.71
%             continue
%         end
%         if mean_I < 0.63 || mean_I > 0.96
%             continue
%         end
%         
        %Filter on avg color.
        % get the avg r,g,b 
        % Filter when every blob r,g,b >=100
        %For tennis ball the color is 198(r),237(g),44(b)
        
%         rgb_mean = zeros(1,3);
%         for k=1:3
%             sum = 0;
%             for j=1:size(rows,2)
%                 sum = sum + curr_img(rows(j) , cols(j) , k);
%             end
%             rgb_mean(k) = sum/size(rows,2);
%         end
        
        %If the candidates color doesn't match.
%         if ~(rgb_mean(1) >= 100 && rgb_mean(2) >= 100 && rgb_mean(3) >= 30)
%             continue;
%         end
        
        %Filter on shape, ellipse fitting.
        % Crop an area around the blob, and fit an ellipse.
        % Add a margin of 2.
        margin = 2;
        min_row = min(rows) - margin;
        min_col = min(cols) - margin;
        
        max_row = max(rows) + margin;
        max_col = max(cols) + margin;
        
        %Make an image of the bound as above.
%         new_image = binary_image .* 255;
%         new_image = insertShape(new_image,'Rectangle',[min_col min_row max_col-min_col max_row-min_row],'LineWidth',1);
%         imshow(new_image);
        %pause(1);
        %cropped_img = binary_image(min_row:max_row , min_col:max_col);
        
        %Fit an ellipse to each of the cropped images:
        new_image = binary_image .* 255;
        % Make sure the bounds are all greater than one, and for 
        if min_row <= 0
            min_row = 1;
        end
        if min_col <=0
            min_col = 1;
        end
        if max_row > s(1)
            max_row = s(1);
        end
        if max_col > s(2)
            max_col = s(2);
        end
        cropped_area = new_image(min_row:max_row , min_col:max_col);
        cropped_area = bwmorph(cropped_area,'remove');
        new_image(min_row:max_row , min_col:max_col) = cropped_area;
        new_image = insertShape(new_image,'Rectangle',[min_col min_row max_col-min_col max_row-min_row],'LineWidth',1);
%        imshow(new_image);
        
        % TO the cropped_area, fit an ellipse
        CC_area = bwconncomp(cropped_area);
        rs = [];
        cls = [];
        for q=1:size(CC_area.PixelIdxList,2)
            vec = CC_area.PixelIdxList(q);
            grpIndx = vec{1};
            [r c] = ind2sub(CC_area.ImageSize ,grpIndx);
            rs = [rs ; r];
            cls = [cls ; c];
        end
        
        if size(cls,1) < 5
            continue;
        end
        
         ellipse_t = fit_ellipse(cls, rs);
        if ~isempty(ellipse_t)
%            disp(ellipse_t.long_axis);
            disp(ellipse_t.short_axis);
            if isempty(ellipse_t.X0_in) || isempty(ellipse_t.Y0)
                continue;
            end
            candidates{count} = Point(min_col + ellipse_t.X0_in, min_row + ellipse_t.Y0_in,k);
            count = count + 1;
            %disp(ellipse_t.long_axis/ellipse_t.short_axis);
        end
        
    end
end