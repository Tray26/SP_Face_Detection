data0 = matfile('feature.mat');

test_output_data = data0.output_mat;
save('temporary.mat', 'test_output_data', '-append');

data1 = matfile('temporary.mat');
append_data = data1.append_data;

new_data = [test_output_data;append_data];
output_mat = new_data;
save('temporary.mat', 'new_data', '-append');
save('feature.mat', 'output_mat', '-append');
%}

data2 = matfile('Copy_of_feature.mat');