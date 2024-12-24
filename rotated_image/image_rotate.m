clear;
clc;
close all;
% 读取图片
% imread 返回一个三维数组，其维度为高度×宽度×颜色通道数（3个颜色通道分别对应RGB
% img = imread('lena.jpg');
img = imread('one_piece.jpg');
imshow(img);
title('original image')

% 旋转角度
degree = 30;

% 求解旋转后图像的尺寸
[m, n, o] = size(img);
% Y = ceil(X) 将 X 的每个元素四舍五入到大于或等于该元素的最接近整数。
new_m = ceil(abs(m*cosd(degree)) + abs(n*sind(degree)));
new_n = ceil(abs(n*cosd(degree)) + abs(m*sind(degree)));

% X = zeros(sz1,...,szN) 返回由零组成的 sz1×...×szN 数组，
% 其中 sz1,...,szN 指示每个维度的大小。例如，zeros(2,3) 将返回一个 2×3 矩阵
new_img_forward = zeros(new_m, new_n, o);
new_img_nnp = zeros(new_m, new_n, 0);
new_img_lp = zeros(new_m, new_n, 0);

%% forward mapping matrices
% 前向映射
m1 = [1 0 0; 0 -1 0; -0.5*n 0.5*m 1];
m2 = [cosd(degree) -sind(degree) 0; sind(degree) cosd(degree) 0; 0 0 1];
m3 = [1 0 0; 0 -1 0; 0.5*new_n 0.5*new_m 1];
m4 = [1 0 0; 0 -1 0; 0.5*n 0.5*m 1];

% 列
for i=1:n
    % 行
    for j=1:m
        new_coordinate = [i j 1]*m1*m2*m3;
        % new_coordinate = [i j 1]*m1*m2*m4;
        % Y = round(X) 将 X 的每个元素四舍五入为最近的整数。在舍入机会均等的情况下，即有元素的十进制小数部分为 0.5（在舍入误差内）时，
        % round 函数会偏离零四舍五入到最接近的具有更大幅值的整数。
        col = max(round(new_coordinate(1)), 1);
        row = max(round(new_coordinate(2)), 1);
        % col = ceil(abs(new_coordinate(1)));
        % row = ceil(abs(new_coordinate(2)));
        % 源图像到已经转换好的图像
        new_img_forward(row, col, 1) = img(j, i, 1);
        new_img_forward(row, col, 2) = img(j, i, 2);
        new_img_forward(row, col, 3) = img(j, i, 3);
    end
end

% figure 使用默认属性值创建一个新的图窗窗口。生成的图窗为当前图窗。
f = figure;
imshow(new_img_forward/255), title('forward mapping')
% 获取该图窗的位置、宽度和高度。
f.Position
%% reverse mapping
% 反向映射
%{
可以采用反向映射的方法，即从旋转后的图像出发，找到对应的原图像的点，
然后将原图像中的灰度值传递过来即可，
这样旋转后的图像的每个像素肯定可以对应到原图像中的一个点，
采取不同策略可以让像素对应地更加准确
%}
rm1 = [1 0 0; 0 -1 0; -0.5*new_n 0.5*new_m 1];
rm2 = [cosd(degree) sind(degree) 0; -sind(degree) cosd(degree) 0; 0 0 1];
rm3 = [1 0 0; 0 -1 0; 0.5*n 0.5*m 1];

% 列 n表示列(i)
for i = 1:new_n
    % 行 m表示行(j)
    for j = 1:new_m
        % 旋转完的坐标系到源图像坐标系
        old_coordinate = [i j 1]*rm1*rm2*rm3;

        % Y = round(X) 将 X 的每个元素四舍五入为最近的整数
        col = round(old_coordinate(1));
        % row：行
        row = round(old_coordinate(2));
        
        non_rounding_col = old_coordinate(1);
        non_rounding_row = old_coordinate(2);
        % 防止边界溢出
        % expr1 || expr2 表示使用逻辑短路行为的逻辑 OR 运算。
        % 如果 expr1 为逻辑值 1 (true)，
        % 将不计算 expr2 的结果。每个表达式的计算结果都必须为标量逻辑值。
        if row < 1 || col < 1 || row > m || col > n
            % j表示行，i表示列
            new_img_nnp(j, i, 1:3) = 0;
            new_img_lp(j, i, 1:3) = 0;
        else
            % 最近邻插值
            new_img_nnp(j ,i, 1) = img(row, col, 1);
            new_img_nnp(j ,i, 2) = img(row, col, 2);
            new_img_nnp(j ,i, 3) = img(row, col, 3);
            
            % 双线性插值（bilinear interpolation
            % % Y = floor(X) 将 X 的每个元素四舍五入到小于或等于该元素的最接近整数
            % left = floor(col);
            % % Y = ceil(X) 将 X 的每个元素四舍五入到大于或等于该元素的最接近整数
            % right = ceil(col);
            % top = floor(row);
            % bottom = ceil(row);
            
            left = max(floor(non_rounding_col), 1);
            right = min(ceil(non_rounding_col), n);
            top = max(floor(non_rounding_row), 1);
            bottom = min(ceil(non_rounding_row), m);
            
            a = non_rounding_col - left;
            b = non_rounding_row - top;

            % a = col - left;
            % b = row -top;

            new_img_lp(j, i, 1) = (1-a)*(1-b)*img(top, left, 1) + a*(1-b)*img(top, right, 1) + ...
                (1-a)*b*img(bottom, left, 1) + a*b*img(bottom, right, 1);
            new_img_lp(j, i, 2) = (1-a)*(1-b)*img(top, left, 2) + a*(1-b)*img(top, right, 2) + ...
                (1-a)*b*img(bottom, left, 2) + a*b*img(bottom, right, 2);
            new_img_lp(j, i, 3) = (1-a)*(1-b)*img(top, left, 3) + a*(1-b)*img(top, right, 3) + ...
                (1-a)*b*img(bottom, left, 3) + a*b*img(bottom, right, 3);
        end
    end
end
figure, imshow(new_img_nnp/255), title('nearest neighbor interpolation');
figure, imshow(new_img_lp/255), title('bilinear interpolation');