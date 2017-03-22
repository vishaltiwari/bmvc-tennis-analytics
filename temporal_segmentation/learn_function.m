run('./vlfeat-0.9.20/toolbox/vl_setup');

m = matfile('segmentdata.mat');

all_C = 0.05:0.05:1;
all_order = [ 2, 3, 4 ];
all_kernel = [ {'KChi2'} ];
max_nice = 0;
best_params = {};

for C=all_C 
    for order = all_order 
        for kernel = all_kernel 
hom.kernel = kernel{:};
hom.order = order;
model_name = [ 'SVM-C', num2str(C), hom.kernel, '_O' , num2str(hom.order), '_model.mat' ];
disp(['params: ', num2str(C), ' ', kernel, ' ', num2str(order) ]);
if exist(model_name, 'file') == 2
    load(model_name);
else
    train_y = m.train_label;
    train_y(train_y == 0) = -1;
    train_set = vl_svmdataset(m.train_data, 'homkermap', hom);
    [w, b, info] = vl_svmtrain(train_set, train_y, C);
    save(model_name, 'w','b','info', '-v7.3');
end

test_y = m.test_label;
test_y(test_y == 0) = -1;
test_set = vl_svmdataset(m.test_data, 'homkermap', hom);
[~,~,~, scores] = vl_svmtrain(test_set, test_y, 0, 'model', w, 'bias', b, 'solver', 'none') ;
scores(scores > 0) = 1;
scores(scores < 0) = -1;

conf_mat = confusionmat(test_y',scores)
if max_nice < conf_mat(2,2)
   max_nice = conf_mat(2,2);
   best_params.C = C;
   best_params.order = order;
   best_params.kernel = kernel;
end
       end
   end
end
best_params