function [v] = worldToImage(K , worldCoords)
    worldCoords = [worldCoords 1];
    imageCoords = K * worldCoords';
    imageCoords(1) = imageCoords(1) / imageCoords(3);
    imageCoords(2) = imageCoords(2) / imageCoords(3);
    v = [imageCoords(1) ; imageCoords(2)];
end