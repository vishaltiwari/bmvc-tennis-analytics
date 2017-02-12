function [cell] = enlargeBlob (binary_image, blob, width , height)
    pixel_loc = blob{1};
    rows = floor((pixel_loc - 1) ./ width) + 1;
    cols = mod((pixel_loc - 1) , width) + 1;
    
    nb_pixels = size(pixel_loc,2);
    
    
    
    %% opening:
    
    
    
    %% manual 4NN anding
    for i=1:nb_pixels
        % get the 8D-NN, and check if its not in the blob.
        for a=-1:1
            for b=-1:1
                row_dash = rows(i) + a;
                col_dash = cols(i) + b;
                
                found_flag = 0;
                if row_dash > 0 && col_dash > 0 && row_dash<=height && cols_dash <= width
                   if binary_image(row_dash,col_dash) == 0 && row_dash + 1 <= width && col_dash + 1 <= height
                        if binary_image(row_dash+1,col_dash+1) && binary_image(row_dash+1,col_dash) && binary_image(row_dash,col_dash+1)
                            
                        end
                   end
                end
                
                % Means that coordinates is not in the blob, but is a
                % neighbour pixel
                
                
                
                
            end
        end
    end
    
    cell{1} = [1;2;3;4];
end