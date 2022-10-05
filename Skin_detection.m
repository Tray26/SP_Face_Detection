test_pic = double(imread('TestImagesForPrograms/14.jpg'));
figure(1)
image(test_pic/255);
%imshow(test_pic/255);

[l,w,c] = size(test_pic);

R=test_pic(:,:,1);
G=test_pic(:,:,2);
B=test_pic(:,:,3);

min_RGB = min(test_pic,[],3);
sum_RGB = sum(test_pic,3);
S = 1-3*min_RGB./sum_RGB;
I = sum_RGB/3;
theta = (R-G+R-B)./(2*sqrt((R-G).^2+(R-B).*(G-B)));
theta = rad2deg(acos(theta));



%max0 = min(min(theta));

mask1 = double(B<=G);
mask2 = double(B>G);

H1 = theta.*mask1;
H2 = (360-theta).*mask2;
H = H1+H2;
max0 = min(min(H));

HSI = zeros(l,w,c);
HSI(:,:,1) = H;
HSI(:,:,2) = S;
HSI(:,:,3) = I;

%直接對HUE做篩選
H0 = abs(H-180);
mask3 = double(H0>120);
test_mask = 255*mask3;
test_mask = hole_filling(test_mask);
mask3 = repmat(test_mask,[1 1 3]);
pic = test_pic.*mask3;
%pic = hsv2rgb(HSI);
%hsv2rgb(HSI);

figure(2)
image(test_mask);
colormap(gray(256));
figure(3)
image(pic/255);
%colormap(gray(256));



function ero = erosion(x)
    ero = min(min(min(min(x, x([1,1:end-1],:)), x([2:end,end],:)), x(:,[1,1:end-1])),x(:,[2:end,end]));
end
function dil = dilation(x)
    dil = max(max(max(max(x, x([1,1:end-1],:)), x([2:end,end],:)), x(:,[1,1:end-1])),x(:,[2:end,end]));
end
function filled_pic = hole_filling(x)
    for i = 1:5
        x = erosion(x);
    end
    for i = 1:5
        x = dilation(x);
    end
    filled_pic = x;
end


