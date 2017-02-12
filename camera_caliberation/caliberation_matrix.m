function [K] = caliberation_matrix(imageCoords , worldCoords)
    
    n = size(worldCoords , 1);
    if ~isequal(size(imageCoords,1), size(worldCoords,1))
        error('Points matrices different sizes');
    end
    if size(worldCoords, 2) ~= 3
        error('Points matrices must have three columns, each for X,Y,Z');
    end
    if n < 6
        error('Need at least 4 matching points');
    end

    M = ones(3,4);
    h = [];
    for i=1:n
        x = worldCoords(i,1);
        y = worldCoords(i,2);
        z = worldCoords(i,3);
        u = imageCoords(i,1);
        v = imageCoords(i,2);
        row1 = [x y z 1 0 0 0 0 -u*x -u*y -u*z -u];
        row2 = [0 0 0 0 x y z 1 -v*x -v*y -v*z -v];
        h = [h ; row1 ; row2];
    end
    if n == 6
        [U, ~, ~] = svd(h');
    else
        [U, ~, ~] = svd(h', 'econ');
    end
    K = (reshape(U(:,12), 4, 3)).';
end