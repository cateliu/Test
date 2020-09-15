% clear
% clc
function [p, Para] = ParameterInMatrix(X,Y,z)
% CCD_x = 320;
% CCD_y = 256;
x = X;
y = Y;
Para_1_1 = sum(sum(x.*x));
Para_1_2 = sum(sum(x.*y));
Para_1_3 = sum(sum(x));
Para_2_1 = Para_1_2;
Para_2_2 = sum(sum(y.*y));
Para_2_3 = sum(sum(y));
Para_3_1 = Para_1_3;
Para_3_2 = Para_2_3;
Para_3_3 = 256*320;
Para = [Para_1_1,Para_1_2, Para_1_3;Para_2_1, Para_2_2, Para_2_3;Para_3_1, Para_3_2, Para_3_3];
Para_in = inv(Para);

b1 = sum(sum(X.*z));
b2 = sum(sum(Y.*z));
b3 = sum(sum(z));
% [b1;b2;b3]*inv(Para)
p = Para\[b1;b2;b3];
% sum(sum(z-(p(1)*x+p(2)*y+p(3)).^2));
% mesh(X,Y,z)
% hold on
% mesh(X,Y,p(1)*x+p(2)*y+p(3))