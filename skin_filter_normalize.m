data = matfile('feature.mat');
output_mat = data.output_mat;

Y_input = output_mat(:,2);
Y_input_mean = mean(Y_input);

Y_input_var = mean((Y_input - Y_input_mean).^2);
Y_input_sd = sqrt(Y_input_var);

Y_output = (Y_input - Y_input_mean)/Y_input_sd;

output_mat(:,2) = Y_output;

training_data = output_mat;

save('feature.mat', 'training_data', '-append');

%data = matfile('feature.mat');
