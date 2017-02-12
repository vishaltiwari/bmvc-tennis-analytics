function [template , templateFeatures , templateDescriptors] = getTemplate(template_in)

% Re-write the function, using a pre-defined set of points where to calculate the sift features and descriptor. 
% First manually mark the points on the image.
% Calculate the features and descriptors at that location of the image.
% return those features and desciptors.

template_img = imread(template_in);
%Use the edges of the template_img
template_img_edge = edge(rgb2gray(template_img),'sobel');   

%Apply hough transform and get the straight lines only.
[H , T, R] = hough(template_img_edge);
P  = houghpeaks(H,20,'threshold',ceil(0.1*max(H(:))));
lines = houghlines(template_img_edge,T,R,P,'FillGap',5,'MinLength',50);
figure ; imshow(template_img_edge);
displayHoughLines(lines);

template = single(template_img_edge);

[templateFeatures , templateDescriptors] = vl_sift(template);

%Display the features.
figure ; imshow(template_img);
displaySIFTFeatures(templateFeatures , templateDescriptors);
end