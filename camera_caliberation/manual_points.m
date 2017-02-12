frame_image = imread('sample_data/camera_caliberation_image.png');
imshow(frame_image)

%[cols_image0 , rows_image0] = ginput();
%[cols_image1 , rows_image1] = ginput();
[cols_image2 , rows_image2] = ginput();
[cols_image3 , rows_image3] = ginput();
[cols_image4 , rows_image4] = ginput();

col_image = (cols_image0 + cols_image1 +cols_image2 +cols_image3 + cols_image4)/5
rows_image = (rows_image0 + rows_image1 +rows_image2 +rows_image3 + rows_image4)/5

