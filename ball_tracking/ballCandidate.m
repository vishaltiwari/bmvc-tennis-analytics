
function [candidates_collection] = ballCandidate (binary_image_boundry , vidHeight,k)
% 1. Get the connected components of the blobs
% 2. Select a one pixel boundery around each blob, and apply enlarging
% 3. Run a sobel filter to get the edges of the blob.
% 4. Run an elipse fitting code on each edges, and get the edges of each
% blob if an elipse was able to get defined.
% 5. For each elipse get the alpha_i, angle between the elipse normal and
% the gradient at that point.
% 6. Filter the elipse and the associated bolb with the following
% parameters: 
%       a. Number of pixels in the blob.
%       b. The alpha value range: [0 0.31] or [0 0.4] or [0 0.5] , less
%       alpha meaning more ball originated blob
%       c. limit of axis length of the elipse (for shape).

% Add the path to the fit_ellipse module:
addpath('/home/vishal/libraries/fit_ellipse')

CC = bwconncomp(binary_image_boundry);
components = CC.PixelIdxList;

new_image = binary_image_boundry .* 255;
count = 1;
%S = regionprops(CC,'Centroid');

ellipse_collection = cell(1,CC.NumObjects);
candidates_collection = cell(1,CC.NumObjects);
for l=1:size(components,2)
    blob_boundry = components(l);
    
    % fit the ellipse in that boundry:
    pixel_loc = blob_boundry{1};
    
    if size(pixel_loc,1) < 5
        continue;
    end
        
    cols = floor((pixel_loc - 1) ./ vidHeight) + 1; % x coordinates
    rows = mod((pixel_loc - 1) , vidHeight) + 1; % y coordinates
    
    ellipse_t = fit_ellipse(cols, rows);
    
    if size(ellipse_t) ~= 0
        if size(ellipse_t.long_axis,1) == 1 && ellipse_t.long_axis <= 20 && ellipse_t.long_axis/ellipse_t.short_axis >=0.5
            %% TODO: Find the alpha value, mean angle difference between the normal and derivate at each pixel point.
            %ellipse_collection{count} = ellipse_t;
            c = Point(ellipse_t.X0_in,ellipse_t.Y0_in,k);
            
            candidates_collection{count} = c;
%             sum = 0;
%             for i=1:size(cols,1)
%                 x0 = cols(i);
%                 y0 = rows(i);
%                 grad = Gdir(y0,x0); %In deg, convert to rad
%                 if grad < 0
%                     grad = 360 + grad;
%                 end
%                 
%                 grad = pi / 180 * grad; % in rad
%                 
%                 alpha_value = getAlphaValue(ellipse_t,x0,y0,grad); % incorrect call, send the gradient of the image as well.
%                 if alpha_value == -1
%                     continue;
%                 end
%                 sum = sum + alpha_value;
%             end
%             
%             mean_alpha = sum / size(cols,1);
%             
%             if mean_alpha <= 1
%                 ellipse_collection{count} = ellipse_t;
%             else
%                 continue;
%             end
            
            %mark the blob:
            %new_image = insertShape(new_image,'circle',[sum(cols)/size(cols,1) , sum(rows)/size(rows,1) , 10] );
            
            if ellipse_t.long_axis <=50
                new_image = insertShape(new_image,'circle',[floor(ellipse_t.X0_in) , floor(ellipse_t.Y0_in) , floor(ellipse_t.long_axis)] );
            else 
                new_image = insertShape(new_image,'circle',[floor(ellipse_t.X0_in) , floor(ellipse_t.Y0_in) , 5] );
            end
            count = count + 1;
        end
    end
    
end
%imshow(new_image);
% Enlarge the connected components:
% enlargeBlob(CC.)
end
