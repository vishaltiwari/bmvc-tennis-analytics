%plot the candidate over time.
X_vec = [];
Y_vec = [];
K_vec = [];
for i=1:size(candidate_per_frame,2)
    candidate_list = candidate_per_frame{i};
    for j=1:size(candidate_list,2)
         if isempty(candidate_list{j})
             continue;
         end
        X_vec = [X_vec candidate_list{j}.x];
        Y_vec = [Y_vec candidate_list{j}.y];
        K_vec = [K_vec candidate_list{j}.k];
    end
end
% Plot the x,y,k in a 3D plot
scatter3(X_vec,Y_vec,K_vec,'x')
