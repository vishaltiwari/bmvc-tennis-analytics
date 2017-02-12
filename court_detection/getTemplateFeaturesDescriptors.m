function [templateFeatures , templateDescriptors , pos] = getTemplateFeaturesDescriptors(template_in)

%Load the manually clicked points.
template_img = imread(template_in);
template_img_edge = edge(imgaussfilt(rgb2gray(template_img),1),'sobel');
[height , width ,dim] = size(template_img);
%load('markedPoints_template_test1.mat');

%Apply hough transform and get the straight lines only.
[H , T, R] = hough(template_img_edge);
P  = houghpeaks(H,20,'threshold',ceil(0.1*max(H(:))));
lines = houghlines(template_img_edge,T,R,P,'FillGap',5,'MinLength',60);
%imshow(template_img);
%displayHoughLines(lines);
court_image = drawLines(lines,width,height,template_img);

%court_image = court_image .* 255;

template_gray = single(imgaussfilt(rgb2gray(court_image),1));

% markedCols, markedRows containing feature locations.
% scales = ones(size(markedCols)) *  10;
% orient = zeros(size(markedCols));
% framesCoords = [markedCols markedRows scales orient];
% 
% [templateFeatures , templateDescriptors] = vl_sift(template_gray , 'frames' , framesCoords' , 'orientations');

[templateFeatures , templateDescriptors] = vl_sift(template_gray,'PeakThresh',10);

%imshow(court_image);
%displaySIFTFeatures(templateFeatures);

%Find the points on the court, and project them.
[cols,rows] = getLinesIntersection(lines,width, height);
pos = [cols rows];
%tempImage = insertShape(template_img,'circle',posM , 'LineWidth',5);

% figure ; imshow(tempImage);

% Or get all the lines, make a court from the line, use that image as
% template.
end