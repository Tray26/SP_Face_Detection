number = '50';
closing_path = sprintf('./closing_result_with_pre_close/%s.png', number);
save_path = sprintf('./skin_filter_final_result/%s.png', number);

pic = double(imread(closing_path));
%{
figure
image(pic)
colormap(gray(256));
%}
labeled_pic = bwlabel(pic);
region_count = max(max(labeled_pic));
%{
figure;
image(labeled_pic)
colormap(gray(256/region_count));
%}
Area = zeros(1,region_count);
for i = 1:region_count
    Area(i) = sum(labeled_pic == i, 'all');
end

[m,n] = size(pic);
pic_area = m*n;
threshold = pic_area/1000;

for i = 1:m
    for j = 1:n
        if pic(i,j) ~= 0
            if Area(labeled_pic(i,j)) <= threshold
                pic(i,j) = 0;
            end
        end
    end
end

figure(1)
image(pic)
colormap(gray(256));
imwrite(pic, save_path)

