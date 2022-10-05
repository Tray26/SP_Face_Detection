clear; clc; clear all;

number = '13';
test_pic = double(imread(sprintf('./TestImagesForPrograms/%s.jpg', number)));
ellipse_parameter = matfile(sprintf('./ellipse_parameters/%s.mat', number)).Z0;
save_map_path = sprintf('./eyemap_indexes/%s.mat', number);

YCC_pic = YCC_transform(test_pic);
eyemapl_result = eyemapl(YCC_pic, 0.05);
eyemapc_result = eyemapc(YCC_pic);
eyemapt_result = eyemapt(YCC_pic, [1;0.2;0.05], 10);


norm_eyemapl = normalization(eyemapl_result);
norm_eyemapc = normalization(eyemapc_result);
norm_eyemapt = normalization(eyemapt_result);

w = [0.45, 0.45, 0.1];
eyemap_ = w(1).*norm_eyemapl + w(2).*norm_eyemapc + w(3)*norm_eyemapt;
h1 = [ellipse_parameter{:,9}];

radius = h1/40;
radius = reshape(radius, [1,1,max(size(radius))]);

radius_mat = repmat(radius, [2*ceil(max(radius))+1, (2*ceil(max(radius))+1), 1]);

eyemap_addi = zeros([size(eyemap_),size(radius_mat,3)]);

[xx,yy] = meshgrid((1:2*ceil(max(radius))+1),(1:2*ceil(max(radius))+1));
circle = (xx-(ceil(max(radius))+1)).^2 + (yy-(ceil(max(radius))+1)).^2;
circle = repmat(circle,[1,1,size(radius,3)]);

cir_matrix = (1./radius_mat).^2.*(circle<=radius_mat.^2);

eyemap_rep = repmat(eyemap_, [1,1,size(radius_mat,3)]);

for i = 1:size(eyemap_addi, 3)
    eyemap_addi(:,:,i) = conv2(eyemap_, cir_matrix(:,:,i), 'same');
end

eyemap0 = zeros(size(eyemap_addi));
eyemap0(:,:,1:end) = (eyemap_addi == max(max(max(max(eyemap_addi, eyemap_addi([1,1:end-1],:,1:end)), eyemap_addi([2:end,end],:, 1:end)), eyemap_addi(:,[1,1:end-1], 1:end)),eyemap_addi(:,[2:end,end], 1:end)));

eyemap1 = ellipse_relation_matching(eyemap0, ellipse_parameter);

eyemap2 = eyemap_addi > 0.2;

eyemap_result = logical(eyemap0 .* eyemap1 .* eyemap2);
for i = 2:size(eyemap_result, 3)
    eyemap_result(:,:,i) = eyemap_result(:,:,i) | eyemap_result(:,:,i-1);
end

sum(eyemap_result(:,:,1), 'all');

% save(save_map_path, "eyemap_result");

[col_axis,row_axis] = find(eyemap_result(:,:,end) == 1);
eyemap_index = [col_axis,row_axis];
save(save_map_path, "eyemap_index");



figure;
image(eyemap2(:,:,2)*255)
colormap(gray(256));

meshgrid(1:size(eyemap0,2),1:size(eyemap0,1));
%}
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

function normalization_result = normalization(matrix)
    normalization_result = (matrix-mean(matrix, 'all'))/std(matrix, 1, 'all');
end

function eyemapt_result = eyemapt(YCC_pic, sigma, L)
    filter_mat = zeros(size(sigma, 1),2*L+1);
    for i = 1:size(sigma, 1)
        filter_mat(i,:) = edge_filter(sigma(i), L);
    end
    
    Y_mat = YCC_pic(:,:,1);
    Y_1x = conv2(Y_mat, filter_mat(1,:), 'same');
    Y_1y = conv2(Y_mat, filter_mat(1,:).', 'same');
    Y_2x = conv2(Y_mat, filter_mat(2,:), 'same');
    Y_2y = conv2(Y_mat, filter_mat(2,:).', 'same');
    Y_3x = conv2(Y_mat, filter_mat(3,:), 'same');
    Y_3y = conv2(Y_mat, filter_mat(3,:).', 'same');
    
    eyemapt_result = max(max(max(max(max(abs(Y_1x),abs(Y_1y)),abs(Y_2x)),abs(Y_2y)),abs(Y_3x)),abs(Y_3y));
end

function filter = edge_filter(sigma, L)
    n = -L:1:L;
    C = 1./sum(exp(-sigma.*[1:L]));
    filter = C.*exp(-sigma.*n).*(n>=1) - C.*exp(-sigma.*abs(n)).*(n<=-1) + 0.*(n==0);
end

function eyemapc_result = eyemapc(YCC_pic)
    cb = YCC_pic(:,:,2);
    cr = YCC_pic(:,:,3);
    cb0 = cb+127.5;
    cr0 = cr+127.5;
    cr1 = max(cr0,[],'all') - cr0;

    eyemapc_result = (cb0.^2+cr1.^2+cb0./cr0)/3;

end

function eyemapl_result = eyemapl(YCC_pic, lambda)
    Y_mat = YCC_pic(:,:,1);
    Y_mat_ero = erosion_n_times(Y_mat, 3);
    Y1 = Y_mat_ero/255;
    eyemapl_result = (1-Y1)./(1+lambda*Y1);
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

function erosion_n = erosion_n_times(x,n)
    erosion_n = x;
    for i = 1:n
        erosion_n = erosion(erosion_n);
    end
end

function ero = erosion(x)
    ero = min(min(min(min(x, x([1,1:end-1],:)), x([2:end,end],:)), x(:,[1,1:end-1])),x(:,[2:end,end]));
end



% function ellipse_parameter = ellipse_matching(pic_order)
%     pic_path = sprintf('./skin_filter_final_result/%s.png', pic_order);
%     pic = double(imread(pic_path));
%     
%     labeled_pic = bwlabel(pic);
%     region_count = max(max(labeled_pic));
% 
% 
%     Area = zeros(1,region_count);
%     for i = 1:region_count
%         Area(i) = sum(labeled_pic == i, 'all');
%     end
%     center_points = zeros(region_count,2);
%     indexes_cell = cell(region_count, 2);
%     Z_Z1 = cell(region_count, 9);
% 
%     for i =1:region_count
%         [m,n] = find(labeled_pic == i);
%         indexes_cell(i, :) = {m,n};
%         m0 = mean(m);
%         n0 = mean(n);
%         center_points(i,:) = [m0,n0];
%         Z_Z1{i,1} = [m-m0,n-n0];
%         Z_Z1{i,2} = [m-m0,n-n0].' * [m-m0,n-n0];
%         [E,V] = eig(Z_Z1{i,2});
%         l1 = V(1,1);
%         l2 = V(2,2);
%         if l1 > l2
%             Z_Z1{i,3} = E(:,1);
%         else
%             Z_Z1{i,3} = E(:,2);
%         end
%         X = [m-m0,n-n0] * E;
%         Z_Z1{i,4} = X;
%         Z_Z1{i,8} = [m-m0,n-n0] * Z_Z1{i,3};
%         m_matrix = zeros(2,2);
%         X1 = X(:,1);
%         X2 = X(:,2);
%         m_matrix(1,1) = mean(abs(X1));
%         m_matrix(1,2) = mean(abs(X2));
%         m_matrix(2,1) = mean(abs(X1).^2);
%         m_matrix(2,2) = mean(abs(X2).^2);
%         Z_Z1{i,5} = m_matrix;
%         
%         ab_matrix = zeros(2,2);
%         ab_matrix(1,:) = 3*pi/4 * m_matrix(1,:);
%         ab_matrix(2,:) = sqrt(4*m_matrix(2,:));
%     
%         a = mean(ab_matrix(:,1));
%         b = mean(ab_matrix(:,2));
%         Z_Z1{i,6} = [a,b];
%     
% %         sum_ = (X1.^2)/a^2 + (X2.^2)/b^2;
%         percentile = sum((X1.^2)/a^2 + (X2.^2)/b^2 <= 1, 'all')/size(X,1);
%         %Z_Z1{i,7} = percentile;
%         ellipse_area = a*b*pi;
%     
%         if percentile >= 0.7 && a/b < 3 && a/b>1/3 && Area(i) > ellipse_area*0.7 && Area(i) < ellipse_area*1.3
%             Z_Z1{i,7} = true;
%         else
%             Z_Z1{i,7} = false;
%             for k = 1:size(n,1)
%                 labeled_pic(m(k), n(k)) = 0;
%                 pic(m(k), n(k)) = 0;
%             end
%         end    
%         Z_Z1{i,9} = max(Z_Z1{i,8}) - min(Z_Z1{i,8});
%     end
%     ellipse_parameter = Z_Z1;
% end

