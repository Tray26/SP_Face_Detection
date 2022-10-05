clear; clc; clear all;

number = '13';
test_pic = double(imread(sprintf('./TestImagesForPrograms/%s.jpg', number)));
ellipse_matching_result = double(imread(sprintf('./ellipse_matching_result/%s.png', number)))/255;
ellipse_parameter = matfile(sprintf('./ellipse_parameters/%s.mat', number)).Z0;
save_mouthmap_path = sprintf('./mouthmap_indexes/%s.mat', number);

%fill the hole on the filter
[m,n] = size(ellipse_matching_result);
k = round(min(m, n)/100)*4;
ellipse_matching_result = closing(ellipse_matching_result, k);
ellipse_matching_result = opening(ellipse_matching_result, k);

YCC_pic = YCC_transform(test_pic);
YCC_pic_filter = YCC_pic.*ellipse_matching_result;
mouthmap_result = mouthmapping(YCC_pic, ellipse_matching_result);

mouthmap_result_norm = normalization(mouthmap_result);
mouthnap_conv_result = convolution(ellipse_parameter, mouthmap_result_norm);
mouthnap_conv_result = normalization(mouthnap_conv_result);

mouthmap0 = zeros(size(mouthnap_conv_result));
mouthmap0(:,:,1:end) = (mouthnap_conv_result ==max(max(max(max(...
    mouthnap_conv_result,...
    mouthnap_conv_result([1,1:end-1],:,1:end)),...
    mouthnap_conv_result([2:end,end],:, 1:end)),...
    mouthnap_conv_result(:,[1,1:end-1], 1:end)),...
    mouthnap_conv_result(:,[2:end,end], 1:end)));

mouthmap1 = ellipse_relation_matching(mouthmap0, ellipse_parameter);
mouthmap2 = (mouthnap_conv_result > 0);

final_mouthmap = logical(mouthmap0 .* mouthmap1 .* mouthmap2);
for i = 2:size(final_mouthmap, 3)
    final_mouthmap(:,:,i) = final_mouthmap(:,:,i) | final_mouthmap(:,:,i-1);
end

final_mouthmap_flat = final_mouthmap(:,:,end);



figure;
image(final_mouthmap_flat*255)
colormap(gray(256));
%}

[col_axis,row_axis] = find(final_mouthmap_flat == 1);
mouthmap_index = [col_axis,row_axis];

save(save_mouthmap_path, "mouthmap_index");


function eyemap1 = ellipse_relation_matching(eyemap0, ellipse_parameter)
    [x_axis,y_axis] = meshgrid(1:size(eyemap0,2),1:size(eyemap0,1), 1:size(eyemap0,3));
    ellipse_center = reshape([ellipse_parameter{:,10}], [2, size([ellipse_parameter{:,10}], 2)/2]);
    ellipse_center = reshape(ellipse_center, [2,1,size(eyemap0,3)]);
    
    x_axis_n = x_axis - ellipse_center(2,:,:);
    y_axis_n = y_axis - ellipse_center(1,:,:);
    
    E = zeros(2,2,size(eyemap0,3));
    V = zeros(2,2,size(eyemap0,3));
    
    x_axis0 = reshape(x_axis_n, [size(x_axis_n, 1)*size(x_axis_n, 2), 1, size(x_axis_n, 3)]);
    y_axis0 = reshape(y_axis_n, [size(y_axis_n, 1)*size(y_axis_n, 2), 1, size(y_axis_n, 3)]);
    
    tran_axis = [x_axis0(:,:,1:end), y_axis0(:,:,1:end)];
    
    for i = 1:size(eyemap0,3)
        [E(:,:,i),V(:,:,i)] = eig(ellipse_parameter{i,2});
        tran_axis(:,:,i) = mtimes(tran_axis(:,:,i), E(:,:,i));
    end
    
    x_axis0 = tran_axis(:,1,:);
    y_axis0 = tran_axis(:,2,:);
    
    x_axis1 = reshape(x_axis0, size(x_axis_n));
    y_axis1 = reshape(y_axis0, size(y_axis_n));
    
    % axis_length = [ellipse_parameter{:,6}];
    axis_length = reshape([ellipse_parameter{:,6}], [2, size([ellipse_parameter{:,6}], 2)/2]);
    axis_length = reshape(axis_length, [2,1,size(tran_axis, 3)]);
    
    a = axis_length(2,:,:);
    b = axis_length(1,:,:);
    
    eyemap1 = ((x_axis1./a).^2 + (y_axis1./b).^2 <= 0.9);
end

function convolution_result = convolution(ellipse_parameter, mouthmap_result_norm)
    height = reshape([ellipse_parameter{:,9}]/25, 1, 1, []);
    width = reshape([ellipse_parameter{:,8}]/4, 1, 1, []);

    height_mat = repmat(height, [2*ceil(max([height; width], [], 'all'))+1, (2*ceil(max([height; width], [], 'all'))+1), 1]);
    width_mat = repmat(width, [2*ceil(max([height; width], [], 'all'))+1, (2*ceil(max([height; width], [], 'all'))+1), 1]);

    mouthmap_addi = repmat(mouthmap_result_norm, 1,1,size(height_mat, 3));

    [xx,yy] = meshgrid((1:2*ceil(max([height; width], [], 'all'))+1),(1:2*ceil(max([height; width], [], 'all'))+1), (1:size(width_mat, 3)));
    xx = xx-ceil(size(xx,1)/2);
    yy = yy-ceil(size(yy,1)/2);
    ellipse = ((yy./height_mat).^2 + (xx./width_mat).^2 <= 1);
%     figure;
%     image(ellipse(:,:,2)/max(ellipse, [], 'all')*255)
%     colormap(gray(256));

    convolution_result = convn(mouthmap_addi, ellipse, 'same');
end

function mouthmap_result = mouthmapping(YCC_pic, ellipse_matching_result)
%     YCC_pic_filter = YCC_pic.*ellipse_matching_result;
    cb = YCC_pic(:,:,2);
    cr = YCC_pic(:,:,3);
    cb0 = cb+127.5;
    cr0 = cr+127.5;
    eta = 0.95*sum((cr0.^2).*ellipse_matching_result, 'all')/ sum((cr0./cb0).*ellipse_matching_result, 'all');

    mouthmap_result = tl(cr0.^2, (cr0.^2-eta*cr0./cb0).^2).*ellipse_matching_result;
end

function tl_result = tl(x,y)
    tl_result = max(x+y-127.5^2, 0);
end

function YCC_pic = YCC_transform(test_pic)
    RGB2YCbCr_tran = [0.299, 0.587, 0.114;...
                       -0.169, -0.331, 0.500;...
                       0.500, -0.419, -0.081];
    test_pic_reshape = reshape(test_pic, [size(test_pic, 1)*size(test_pic, 2), 3]).';
    
    YCC_pic_reshape = zeros(size(test_pic_reshape));
    YCC_pic_reshape(:,1:end) = RGB2YCbCr_tran * test_pic_reshape(:,1:end);
    YCC_pic_reshape = YCC_pic_reshape.';
    
    YCC_pic = reshape(YCC_pic_reshape, size(test_pic));
end

function normalization_result = normalization(matrix)
    normalization_result = (matrix-mean(matrix, 'all'))/std(matrix, 1, 'all');
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


% eyemap_result = matfile(sprintf('./eyemap_results/%s.mat', number)).eyemap_result;
% eye_map_path = sprintf('./eyemap_results/%s.mat', number);
% figure;
% image(mouthmap_result_norm/max(mouthmap_result_norm, [], 'all')*255)
% colormap(gray(256));
% YCC_pic_filter(:,:,2);
% cb = YCC_pic(:,:,2);
% cr = YCC_pic(:,:,3);
% cb0 = cb+127.5;
% cr0 = cr+127.5;
% cr0_filtered = cr0.*ellipse_matching_result;
% cb0_filtered = cb0.*ellipse_matching_result;
% cr0_filtered./cb0_filtered;
% 
% eta = sum((cr0.^2).*ellipse_matching_result, 'all')/ sum((cr0./cb0).*ellipse_matching_result, 'all');
% figure;
% image(ellipse_matching_result*255)
% colormap(gray(256));
%{
height = reshape([ellipse_parameter{:,9}]/25, 1, 1, []);
width = reshape([ellipse_parameter{:,8}]/4, 1, 1, []);

height_mat = repmat(height, [2*ceil(max([height; width], [], 'all'))+1, (2*ceil(max([height; width], [], 'all'))+1), 1]);
width_mat = repmat(width, [2*ceil(max([height; width], [], 'all'))+1, (2*ceil(max([height; width], [], 'all'))+1), 1]);

mouthmap_addi = zeros([size(mouthmap_result_norm),size(height_mat,3)]);

[xx,yy] = meshgrid((1:2*ceil(max([height; width], [], 'all'))+1),(1:2*ceil(max([height; width], [], 'all'))+1), (1:size(width_mat, 3)));
xx = xx-ceil(size(xx,1)/2);
yy = yy-ceil(size(yy,1)/2);
ellipse = ((yy./height_mat).^2 + (xx./width_mat).^2 <= 1);
convolution_result = convn(mouthmap_addi, ellipse, 'same');
%}
