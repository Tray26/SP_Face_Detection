clear; clc; clear all;

number = '01';
test_pic = double(imread(sprintf('./TestImagesForPrograms/%s.jpg', number)));
eyemap_index = matfile(sprintf('./eyemap_indexes/%s.mat', number)).eyemap_index;
mouthmap_index = matfile(sprintf('./mouthmap_indexes/%s.mat', number)).mouthmap_index;
ellipse_parameter = matfile(sprintf('./ellipse_parameters/%s.mat', number)).Z0;

ellipses = ellipse_finding(test_pic, ellipse_parameter);
eyes_mouth_pair = eyes_mouth_pairing(eyemap_index, mouthmap_index);
[ellipse_index, eyes_mouth_pair] = ellip_filter(eyes_mouth_pair, ellipses);
face_dectect_result = feature_mapping(eyes_mouth_pair, ellipse_index, ellipse_parameter);

face_dectection_result = (face_dectect_result ~= 0);
%}





function feature_map_3 = principle_axis_matching(vector_cd, ellipse_parameter, ellipse_index)
    principal_axises = [ellipse_parameter{:,3}].';
    principal_axis = principal_axises(ellipse_index,:);

    Dot = dot(principal_axis, vector_cd, 2);
    norm_principal_axis = sqrt(principal_axis(:,1).^2+principal_axis(:,2).^2);
    norm_cd = sqrt(vector_cd(:,1).^2+vector_cd(:,2).^2);
    angle = acos(Dot./(norm_principal_axis.*norm_cd))*360/(2*pi);
    feature_map_3 = (angle<5 | angle>175);
end

function eyes_mouth_pair = eyes_mouth_pairing(eyemap_index, mouthmap_index)
    eyes_mouth_pair = [];
    for i =1:size(eyemap_index, 1)-1
        e2 = eyemap_index(i+1:end, :);
        e1 = repmat(eyemap_index(i,:), size(e2, 1), 1);
        eyes = cat(3, e1, e2);
        eyes_rep = repmat(eyes, size(mouthmap_index, 1), 1, 1);
        mouthmap_index_rep = repmat(mouthmap_index.', size(eyes, 1), 1, 1);
        mouthmap_index_rep = reshape(mouthmap_index_rep, 2, []).';
        new_eyes_mouth_pair = cat(3, eyes_rep, mouthmap_index_rep);
        eyes_mouth_pair = [eyes_mouth_pair;new_eyes_mouth_pair];
    end
end

function [ellipse_index, ellip_filtered_eyes_mouth_pair] = ellip_filter(eyes_mouth_pair, ellipses)
    ellipse_index = [];
    ellip_filtered_eyes_mouth_pair = [];
    for i = 1:size(eyes_mouth_pair, 1)
        pair = reshape(eyes_mouth_pair(i,:,:), 2,[]).';
        ellip_index = find(ellipses(pair(1,1),pair(1,2),:) == 1 & ...
                                  ellipses(pair(2,1), pair(2,2), :) == 1 & ...
                                  ellipses(pair(3,1), pair(3,2), :) == 1);
        if(size(ellip_index)~=0)
            ellipse_index(end+1,1) = ellip_index;
            ellip_filtered_eyes_mouth_pair(end+1,:,:) = eyes_mouth_pair(i,:,:);
        end
    end
end

function ellipses = ellipse_finding(input_pic, ellipse_parameter)
    [x_axis,y_axis] = meshgrid(1:size(input_pic,2),1:size(input_pic,1), 1:size(ellipse_parameter,1));
    ellipse_center = reshape([ellipse_parameter{:,10}], [2, size([ellipse_parameter{:,10}], 2)/2]);
    ellipse_center = reshape(ellipse_center, [2,1,size(ellipse_parameter,1)]);
    
    x_axis_n = x_axis - ellipse_center(2,:,:);
    y_axis_n = y_axis - ellipse_center(1,:,:);
    
    E = zeros(2,2,size(ellipse_parameter,1));
    V = zeros(2,2,size(ellipse_parameter,1));
    
    x_axis0 = reshape(x_axis_n, [size(x_axis_n, 1)*size(x_axis_n, 2), 1, size(x_axis_n, 3)]);
    y_axis0 = reshape(y_axis_n, [size(y_axis_n, 1)*size(y_axis_n, 2), 1, size(y_axis_n, 3)]);
    
    tran_axis = [x_axis0(:,:,1:end), y_axis0(:,:,1:end)];
    
    for i = 1:size(ellipse_parameter,1)
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
    
    ellipses = ((x_axis1./a).^2 + (y_axis1./b).^2 <= 0.9);

end

function face_dectect_result = feature_mapping(eyes_mouth_pair, ellipse_index, ellipse_parameter)
    D_mid_eye = mean(eyes_mouth_pair(1:end,:, 1:2), 3);
    vector_ab = eyes_mouth_pair(:, :, 1) - eyes_mouth_pair(:, :, 2);
    vector_cd = D_mid_eye - eyes_mouth_pair(:,:,3);

    feature1 = (eyes_mouth_pair(1:end, 1, 1) > eyes_mouth_pair(1:end, 1, 3) &...
                     eyes_mouth_pair(1:end, 1, 2) > eyes_mouth_pair(1:end, 1, 3));
    feature2 = feature_mapping_step2(vector_ab, vector_cd);
    feature3 = principle_axis_matching(vector_cd, ellipse_parameter, ellipse_index);
    feature4 = feature_mapping_step4(vector_ab, vector_cd);

    face_dectect_result = sum(feature1 & feature2 & feature3 & feature4);
end

function feature_map_4 = feature_mapping_step4(vector_ab, vector_cd)
    norm_ab = sqrt(vector_ab(:,1).^2+vector_ab(:,2).^2);
    norm_cd = sqrt(vector_cd(:,1).^2+vector_cd(:,2).^2);
    distance_ratio = norm_ab./norm_cd;
    feature_map_4 = (distance_ratio < 2.5 & distance_ratio > 0.4);
end

function feature_map_2 = feature_mapping_step2(vector_ab, vector_cd)
    Dot = dot(vector_ab, vector_cd, 2);
    norm_ab = sqrt(vector_ab(:,1).^2+vector_ab(:,2).^2);
    norm_cd = sqrt(vector_cd(:,1).^2+vector_cd(:,2).^2);
    angle = acos(Dot./(norm_ab.*norm_cd))*360/(2*pi);
    feature_map_2 = (angle<95 & angle>85);
end



%{
dectect = zeros(size(eyes_mouth_pair, 1), 4);

dectect(:, 1) = (eyes_mouth_pair(1:end, 1, 1) > eyes_mouth_pair(1:end, 1, 3) &...
                 eyes_mouth_pair(1:end, 1, 2) > eyes_mouth_pair(1:end, 1, 3));

D_mid_eye = mean(eyes_mouth_pair(1:end,:, 1:2), 3);

vector_ab = eyes_mouth_pair(:, :, 1) - eyes_mouth_pair(:, :, 2);
vector_cd = D_mid_eye - eyes_mouth_pair(:,:,3);

feature_map_2 = feature_mapping_step2(vector_ab, vector_cd);
dectect(:,2) = feature_map_2;

feature_map_3 = principle_axis_matching(vector_cd, ellipse_parameter, ellipse_index);
dectect(:,3) = feature_map_3;

dectect(:,4) = feature_mapping_step4(vector_ab, vector_cd);



principal_axises = [ellipse_parameter{:,3}].';
principal_axis = principal_axises(ellipse_index,:);

Dot = dot(principal_axis, vector_cd, 2);

norm_principal_axis = sqrt(principal_axis(:,1).^2+principal_axis(:,2).^2);
norm_cd = sqrt(vector_cd(:,1).^2+vector_cd(:,2).^2);
norm_ab = sqrt(vector_ab(:,1).^2+vector_ab(:,2).^2);
angle = acos(Dot./(norm_principal_axis.*norm_cd))*360/(2*pi);
feature_map_3 = (angle<5 | angle>175);

% ellipse_index = [];
% ellip_filtered_eyes_mouth_pair = [];
% 
% for i = 1:size(eyes_mouth_pair, 1)
%     pair = reshape(eyes_mouth_pair(i,:,:), 2,[]).';
%     ellip_index = find(ellipses(pair(1,1),pair(1,2),:) == 1 & ...
%                               ellipses(pair(2,1), pair(2,2), :) == 1 & ...
%                               ellipses(pair(3,1), pair(3,2), :) == 1);
%     if(size(ellip_index)~=0)
%         ellipse_index(end+1,1) = ellip_index;
%         ellip_filtered_eyes_mouth_pair(end+1,:,:) = eyes_mouth_pair(i,:,:);
%     end
% end
% eyes_mouth_pair_fil = ellip_filtered_eyes_mouth_pair;


for i = 1:size(ellipses, 3)
%     t = find(ellipses(pair(1,1),pair(1,2),i) == 1);
    for j = 1:size(eyes_mouth_pair, 1)

    end
end

% pair = reshape(eyes_mouth_pair(1,:,:), [],2);
% index = [pair(1,:), 1];
% t = ellipses(pair(1,1),pair(1,2),2);

for i = 1:size(eyes_mouth_pair, 1)
    pair = reshape(eyes_mouth_pair(i,:,:), [],2);
    t = ellipses([pair(1,:)]);
    disp(t);
    for j = 1:size(ellipses, 3)

    end
end


    for j = i+1:size(eyemap_index, 1)
        for k =1:size(mouthmap_index, 1)
            eye1 = eyemap_index(i,:);
            eye2 = eyemap_index(j,:);
            mouth = mouthmap_index(k,:);
        end
    end

% eye1_col = eyes_mouth_pair(1:end, 1, 1);
% eye2_col = eyes_mouth_pair(1:end, 1, 2);
% mouth_col = eyes_mouth_pair(1:end, 1, 3);
%}