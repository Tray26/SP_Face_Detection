load_data = matfile('feature.mat');
model = load_data.model;

number = "50";
input_path = sprintf('./TestImagesForPrograms/%s.jpg', number);
filtered_path = sprintf('./skin_filter_results/%s.png', number);
opening_path = sprintf('./opening_result_without_pre_close/%s.png', number);
closing_path = sprintf('./closing_result_without_pre_close/%s.png', number);

test_pic = double(imread(input_path));

test_data = pic_to_test(test_pic);
normalized_test_data = Y_normalize(test_data);

test_feature = sparse(normalized_test_data);
test_label = ones(size(test_feature, 1), 1);

[predicted_label, accuracy, decision_values] = svmpredict(test_label, test_feature, model);

[w,l,c] = size(test_pic);
filtered_pic = reshape(predicted_label, [w,l]);
filtered_pic = filtered_pic * 255;

figure(1)
image(filtered_pic)
colormap(gray(256));
%imwrite(filtered_pic, filtered_path)



correction_points = 2;
t = ginput(correction_points);
t_round = round(t);

%label = ones(correction_points, 1);
label = zeros(correction_points, 1);
%label = [0;1];

RGB_matrix = zeros(correction_points, 3);

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
%}


pre_data = matfile('temporary.mat');
append_data = pre_data.append_data;

append_data(end+1:end+correction_points, :) = [label, YCbCr_matrix];
save('temporary.mat', 'append_data', '-append');
%}
%{
append_data = [label, YCbCr_matrix];
save('temporary.mat', 'append_data');
%}


[m,n] = size(filtered_pic);
k = round(min(m, n)/100);

filtered_pic = closing(filtered_pic, 2);


opened_pic = opening(filtered_pic, k);
figure(2)
image(opened_pic)
colormap(gray(256));
%imwrite(opened_pic, opening_path)

closed_pic = closing(opened_pic, k);
figure(3)
image(closed_pic)
colormap(gray(256));
%imwrite(closed_pic, closing_path)
%}


function YCbCr_matrix = pic_to_test(test_pic)
    [w,l,c] = size(test_pic);
    RGB = reshape(test_pic, [w*l, 1, c]);
    RGB_mat = zeros(size(RGB, 1), c);
    RGB_mat(1:end, :) = RGB(1:end,1, :);
    RGB_mat = RGB_mat .';
    
    RGB_tran_matrix = [0.299, 0.587, 0.114;...
                       -0.169, -0.331, 0.500;...
                       0.500, -0.419, -0.081];

    YCbCr_matrix = zeros(size(RGB_mat));
    YCbCr_matrix(:,1:end) = mtimes(RGB_tran_matrix,RGB_mat(:,1:end));
    
    YCbCr_matrix = YCbCr_matrix .';
end

function normalized_data = Y_normalize(test_data)
    Y_input = test_data(:,1);
    Y_input_mean = mean(Y_input);
    Y_input_var = mean((Y_input - Y_input_mean).^2);
    Y_input_sd = sqrt(Y_input_var);
    Y_output = (Y_input - Y_input_mean)/Y_input_sd;
    normalized_data = test_data;
    normalized_data(:,1) = Y_output;
end

function ero = erosion(x)
    ero = min(min(min(min(x, x([1,1:end-1],:)), x([2:end,end],:)), x(:,[1,1:end-1])),x(:,[2:end,end]));
end

function dil = dilation(x)
    dil = max(max(max(max(x, x([1,1:end-1],:)), x([2:end,end],:)), x(:,[1,1:end-1])),x(:,[2:end,end]));
end

function opening_pic = opening(pic, k)
    for i = 1:k
        pic = erosion(pic);
    end
    for i = 1:k
        pic = dilation(pic);
    end
    opening_pic = pic;
end

function closing_pic = closing(pic, k)
    for i = 1:k
        pic = dilation(pic);
    end
    for i = 1:k
        pic = erosion(pic);
    end
    closing_pic = pic;
end
