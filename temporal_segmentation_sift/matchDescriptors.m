function [flag , tform] = matchDescriptors(templateFatures, templateDescriptor, frameFeatures, frameDescriptors , frame_image , tem_in)
[matches, scores] = vl_ubcmatch(templateDescriptor, frameDescriptors, 3);
flag = 1;
disp(size(matches,2));
if size(matches,2) <= 4 %Changes from == 0 to <= 8, from matching frames in australian open, ~60 matchess
    flag = 0;
    tform=[];
    return;
else %remove else if I want to globally do a motion compensation.
    flag = 1;
    tform = [];
    return;
end
% See the matching points. This is not being 
f_template = templateFatures(:,matches(1,:));
f_frame = frameFeatures(:,matches(2,:));

% Display matching points
%temp = imread(tem_in);
%figure ; imshow(temp);
%displaySIFTFeatures(f_template);
%figure ; imshow(frame_image);
%displaySIFTFeatures(f_frame);

template_coords = f_template(1:2,:);
frame_coords = f_frame(1:2,:);

[tform , inlierpoints1 , inlierpoints2 , status] = ... 
    estimateGeometricTransform(template_coords' , frame_coords' , 'affine');


if status == 1 || status == 2
    disp('Not enough matching points or Not enough inliers points');
    flag = 0;
    return;
end

% if size(inlierpoints1,2) <= 20
%     flag = 0;
%     return;
% end

%% Display matching inlines. FOR DEBUGGING
% inliner1Circles = [inlierpoints1 5*ones(size(inlierpoints1,1),1)];
% inliner2Circles = [inlierpoints2 5*ones(size(inlierpoints1,1),1)];

% temp = imread(tem_in);
% temp = insertShape(temp, 'circle',inliner1Circles);
% imshow(temp);
% 
% frame_image2 = insertShape(frame_image , 'circle' , inliner2Circles);
% imshow(frame_image2);

end