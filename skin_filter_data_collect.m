test_pic = double(imread('./TestImagesForPrograms/5.jpg'));
image(test_pic/255);
sample_points = 2;
t = ginput(sample_points);
t_round = round(t);

label = [0;1];

s = size(t_round);

RGB_matrix = zeros(sample_points, 3);

for i = 1:size(t_round, 1)
    RGB_matrix(i,:) = test_pic(t_round(i, 2), t_round(i, 1), :);
end

RGB_matrix = RGB_matrix.';

RGB_tran_matrix = [0.299, 0.587, 0.114;...
                   -0.169, -0.331, 0.500;...
                   0.500, -0.419, -0.081];

YCbCr_matrix = zeros(size(RGB_matrix));
YCbCr_matrix(:,1:end) = mtimes(RGB_tran_matrix,RGB_matrix(:,1:end));
YCbCr_matrix = YCbCr_matrix.';


example = matfile('feature.mat');
output_mat = example.output_mat;

output_mat(end+1:end+2,:) = [label, YCbCr_matrix];
%}
%output_mat = [label, YCbCr_matrix];

save('feature.mat', "output_mat");
%}
