clear; clc; clear all;
number = '06';
pic_path = sprintf('./skin_filter_final_result/%s.png', number);
save_path = sprintf('./ellipse_matching_result/%s.png', number);
parameter_path = sprintf('./ellipse_parameters/%s.mat', number);

pic = double(imread(pic_path));
% [m,n] = size(pic);
% pic = opening(pic, round(min(m, n)));

labeled_pic = bwlabel(pic);
region_count = max(max(labeled_pic));

center_points = zeros(region_count,2);
indexes_cell = cell(region_count, 2);
Z_Z1 = cell(region_count, 10);

Area = zeros(1,region_count);
for i = 1:region_count
    Area(i) = sum(labeled_pic == i, 'all');
end

for i =1:region_count
    [m,n] = find(labeled_pic == i);
    indexes_cell(i, :) = {m,n};
    m0 = mean(m);
    n0 = mean(n);
    center_points(i,:) = [m0,n0];
    Z_Z1{i,1} = [m-m0,n-n0];                            %Z_Z1{i,1} = All points in the region
    Z_Z1{i,2} = [m-m0,n-n0].' * [m-m0,n-n0];            %Z_Z1{i,2} = Z1
    Z_Z1{i,10} = [m0,n0];
    [E,V] = eig(Z_Z1{i,2});
    l1 = V(1,1);
    l2 = V(2,2);
    if l1 > l2                                          %Z_Z1{i,3} = principal axis
        Z_Z1{i,3} = E(:,1);
        short_axis = E(:,2);
        short_index = [m-m0,n-n0] * short_axis;
        Z_Z1{i,8} = max(short_index) - min(short_index);
    else
        Z_Z1{i,3} = E(:,2);
        short_axis = E(:,1);
        short_index = [m-m0,n-n0] * short_axis;
        Z_Z1{i,8} = max(short_index) - min(short_index);
    end
    X = [m-m0,n-n0] * E;
    Z_Z1{i,4} = X;                                      %Z_Z1{i,4} = location after axis-trans.
%     Z_Z1{i,8} = [m-m0,n-n0] * Z_Z1{i,3};
    m_matrix = zeros(2,2);
    X1 = X(:,1);
    X2 = X(:,2);
    m_matrix(1,1) = mean(abs(X1));
    m_matrix(1,2) = mean(abs(X2));
    m_matrix(2,1) = mean(abs(X1).^2);
    m_matrix(2,2) = mean(abs(X2).^2);
    Z_Z1{i,5} = m_matrix;
    
    ab_matrix = zeros(2,2);
    ab_matrix(1,:) = 3*pi/4 * m_matrix(1,:);
    ab_matrix(2,:) = sqrt(4*m_matrix(2,:));

    a = mean(ab_matrix(:,1));
    b = mean(ab_matrix(:,2));
    Z_Z1{i,6} = [a,b];

    percentile = sum((X1.^2)/a^2 + (X2.^2)/b^2 <= 1, 'all')/size(X,1);
    %Z_Z1{i,7} = percentile;
    ellipse_area = a*b*pi;

    if percentile >= 0.7 && a/b < 3 && a/b>1/3 && Area(i) > ellipse_area*0.7 && Area(i) < ellipse_area*1.3
        Z_Z1{i,7} = true;
    else
        Z_Z1{i,7} = false;
        for k = 1:size(n,1)
            labeled_pic(m(k), n(k)) = 0;
            pic(m(k), n(k)) = 0;
        end
    end
    Z_Z1{i,9} = max([m-m0,n-n0] * Z_Z1{i,3}) - min([m-m0,n-n0] * Z_Z1{i,3});

    
end

Z0 = {};
for i = 1:size(Z_Z1,1)
    if Z_Z1{i,7} == 1
        Z0(end+1,:) = Z_Z1(i,:);
    end
end
% Z0 = Z0(2:end,:);


save(parameter_path, 'Z0');

figure(1)
image(pic)
colormap(gray(256));
imwrite(pic, save_path);
%}
