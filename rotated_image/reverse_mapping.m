% reverse mapping matrices
% 其中，m表示行(j)，n表示列(i)
%{
rm1:仿射变换矩阵，用于处理图像的平移和镜像翻转（关于y轴）
第一行[1 0 0]保持x坐标不变
第二行[0 -1 0]实现y坐标的镜像反转
第三行[-0.5*new_n 0.5*new_m 1]进行图像平移

rm2:是一个旋转矩阵，旋转角度为 degree，单位是度，因为cosd(degree)用于计算角度的余弦和正弦值
%}
rm1 = [1 0 0; 0 -1 0; -0.5*new_n 0.5*new_m 1];
rm2 = [cosd(degree) sind(degree) 0; -sind(degree) cosd(degree) 0; 0 0 1];
rm3 = [1 0 0; 0 -1 0; 0.5*n 0.5*m 1];

% 代表x方向
for i = 1:new_n
    % 代表y方向
    for j = 1:new_m
        % 旋转后的图像坐标系转到无旋转的图像坐标系，即将转换完的图像中的坐标一一对应到原来的源图像坐标系，以便进行颜色和各种信息的匹配
        % 遍历新图像中的每一个像素点，i为列索引（水平位置），j为行索引（垂直位置）
        old_coordinate = [i j 1] * rm1 * rm2 * rm3;
        % old_coordinate(1)：这表示 old_coordinate 向量的第一个元素，即变换后的 x 坐标
        % col = round(old_coordinate(1))：这行代码的作用是将变换后的 x 坐标（old_coordinate(1)）四舍五入到最接近的整数，
        % 并将结果赋值给变量 col。在图像处理中，col 通常表示像素的列索引（即水平位置）。
        % 在舍入机会均等的情况下，即有元素的十进制小数部分为 0.5（在舍入误差内）时，
        % round 函数会偏离零四舍五入到最接近的具有更大幅值的整数。
        col = round(old_coordinate(1));
        row = round(old_coordinate(2));
        % 防止超出边界
        % 检查row和col是否超出原图的边界（即行列数m和n),如果超出边界，则将变换后图像中的当前像素位置（new_img_nnp(j, i) 和 new_img_lp(j, i)）赋值为 0，表示背景或透明
        if row < 1 || col <1 || row>m || col > n
            new_img_nnp(j, i) = 0;
            new_img_lp(j, i) = 0;
        else
            
        end

    end

end

