function [templateFeatures , templateDescriptors , pos] = getTemplateFeaturesDescriptors(template_in)

%NOTE:For each video, need to extract a template of the court.
%using template as it is.

%  binSize = 8 ;
%  magnif = 3 ;
template_img = imread(template_in);
template_img = imresize(template_img , 0.10);
template_gray = single(imgaussfilt(rgb2gray(template_img),1));

% [templateFeatures , templateDescriptors] = vl_dsift(template_gray , 'size', binSize); % For synthetic court/grass court templates(aus/us open)
% templateFeatures(3,:) = binSize/magnif ;
% templateFeatures(4,:) = 0 ;
[templateFeatures , templateDescriptors] = vl_sift(template_gray , 'PeakThresh',0); % For synthetic court/grass court templates(aus/us open)
%[templateFeatures , templateDescriptors] = vl_sift(template_gray , 'PeakThresh',1); % For clay court templates

%[height , width] = size(template_img);

% Find the intersecting lines
%template_img_edge = edge(imgaussfilt(rgb2gray(template_img),1));
%[H , T, R] = hough(template_img_edge);
%P  = houghpeaks(H,20,'threshold',ceil(0.1*max(H(:))));
%lines = houghlines(template_img_edge,T,R,P,'FillGap',5,'MinLength',60);
%[cols,rows] = getLinesIntersection(lines,width, height);
%pos = [cols rows];
pos=[];

figure; imshow(template_img);
displaySIFTFeatures(templateFeatures);

end