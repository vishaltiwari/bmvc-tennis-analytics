function [frameFeatures , frameDescriptors] = getFrameFeatureDescriptor(frame_image)

[height, width, dim] = size(frame_image);
frame_image_edge = edge(imgaussfilt(rgb2gray(frame_image),1),'sobel');

court_image = frame_image;

%Apply hough transform and get the straight lines only.
% [H , T, R] = hough(frame_image_edge);
% P  = houghpeaks(H,20,'threshold',ceil(0.1*max(H(:))));
% lines = houghlines(frame_image_edge,T,R,P,'FillGap',5,'MinLength',60);
%imshow(frame_image);
%displayHoughLines(lines);
% court_image = drawLines(lines,width,height,frame_image);
%court_image = court_image .* 255;
% imshow(frame_image);
% hold on;
% displayHoughLines(lines);

% [cols,rows] = getLinesIntersection(lines,width, height);
% posM = [cols rows 3*ones(size(cols))];
% tempImage = insertShape(frame_image,'circle',posM , 'LineWidth',5);
% imshow(tempImage);
% Find a) interestions, and find b) features, descriptors.

frame_image_gray = single(imgaussfilt(rgb2gray(court_image),1));
% scales = ones(size(cols)) *  7;
% orient = zeros(size(cols));
% framesCoords = [cols rows scales orient];
% framesCoords = unique(framesCoords,'rows');

% [frameFeatures , frameDescriptors] = vl_sift(frame_image_gray , 'frames' , framesCoords' , 'orientations');
[frameFeatures , frameDescriptors] = vl_sift(frame_image_gray,'PeakThresh',0.5);
%figure ; imshow(court_image);
%displaySIFTFeatures(frameFeatures);
end