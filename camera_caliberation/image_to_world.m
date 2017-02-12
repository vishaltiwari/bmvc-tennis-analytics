img = imread('sample_data/camera_caliberation_image.png');
[tform , inlierpoints1 , inlierpoints2 , status] = ... 
    estimateGeometricTransform(imageCoordinatePoints2 , worldCoordinatePoints , 'projective');

%outputView = imref2d(size(img));
outputImage = imwarp(img,tform,'cubic');
imshow(img);
figure ; imshow(outputImage);
