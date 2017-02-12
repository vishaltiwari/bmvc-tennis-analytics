function [v] = getImagePoints(XYZ_sim)

    load('../camera_caliberation/Camera_matrix.mat');

    % use the K matrix to get the image coordinates:
    n = size(XYZ_sim,1);
    worldCoords = [XYZ_sim(:,1:3) ones(n,1)];
    imageCoords = K * worldCoords';
    imageCoords(1,:) = imageCoords(1,:) ./ imageCoords(3,:);
    imageCoords(2,:) = imageCoords(2,:) ./ imageCoords(3,:);
    v = [imageCoords(1,:)'  imageCoords(2,:)'];

end