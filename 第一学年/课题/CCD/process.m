clear
load('Report2_4.mat')

x = 1:320;
y = 1:256;
[X,Y] = meshgrid(x-160, y-128);
iavg = (I(:,:,20)+I(:,:,40))/2;
ag = angle(I(:,:,40)-iavg);
ag1 = ag';
mesh(X,Y,-ag1)
ylim([-128 128])
xlim([-160 160])
yticks([-128,0,128])
xticks([-160,0,160])
xlabel('x轴像素')
ylabel('y轴像素')
title('CCD上相位分布')
X1 = X*30e-6;
Y1 = Y*30e-6;
p = ParameterInMatrix(X1,Y1,ag1)
lambda = 1064e-9;
k = 2*pi/lambda;

p(1)/30e-6/k
function [p, Para] = ParameterInMatrix(x,y,z)
% CCD_x = 320;
% CCD_y = 256;
[a b] = size(x);
Para_1_1 = sum(sum(x.*x));
Para_1_2 = sum(sum(x.*y));
Para_1_3 = sum(sum(x));
Para_2_1 = Para_1_2;
Para_2_2 = sum(sum(y.*y));
Para_2_3 = sum(sum(y));
Para_3_1 = Para_1_3;
Para_3_2 = Para_2_3;
Para_3_3 = a*b;
Para = [Para_1_1,Para_1_2, Para_1_3;Para_2_1, Para_2_2, Para_2_3;Para_3_1, Para_3_2, Para_3_3];


b1 = sum(sum(x.*z));
b2 = sum(sum(y.*z));
b3 = sum(sum(z));
% [b1;b2;b3]*inv(Para)
p = Para\[b1;b2;b3];

sqrt(norm(z-p(1)*x-p(2)*y-p(3)))
end
