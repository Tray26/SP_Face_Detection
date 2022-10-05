data = matfile('feature.mat');
training_data = data.training_data;

labels = training_data(:,1);
features = training_data(:,2:end);
fs = sparse(features);
libsvmwrite('svm', labels, fs);
[y,x] = libsvmread('svm');

train_label = y;
train_data = x;

model = svmtrain(train_label, train_data);

save("feature.mat", 'model', '-append')
