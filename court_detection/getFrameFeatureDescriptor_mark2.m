function [frameFeatures , frameDescriptors] = getFrameFeatureDescriptor_mark2(frame_image)
frame_gray = single(imgaussfilt(rgb2gray(frame_image),1));
[frameFeatures , frameDescriptors] = vl_sift(frame_gray,'PeakThresh',5); % for syntheic/grass courts
%[frameFeatures , frameDescriptors] = vl_sift(frame_gray,'PeakThresh',1); % clay courts. Get more features.
% imshow(frame_image);
% displaySIFTFeatures(frameFeatures);
end