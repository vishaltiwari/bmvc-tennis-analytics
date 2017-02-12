function [] = SwapPaths(paths,weights,lengths,i,j)
    %Swap
    tmp_p = paths(i);
    tmp_w = weights(i);
    tmp_l = lengths(i);
    paths(i) = paths(j);
    weights(i) = weights(j);
    lengths(i) = lengths(j);
    paths(j) = tmp_p;
    weights(j) = tmp_w;
    lengths(j) = tmp_l;
end