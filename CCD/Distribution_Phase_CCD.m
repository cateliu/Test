%% CCD的计算
%干涉的角度确定，绕y轴转delta。alpha = 0；delta = 0.2mrad；
%CCD的微元面积为30e-6 X 30e-6；阵列为 320x256；
% 两道高斯光入射，入射参数相同，仅角度不同。
clear
clc
format long
delta = 2e-3;

w_0 = 2e-3;
lambda = 1064e-9;
k = 2*pi/lambda;
z_r = pi*1*w_0^2/lambda;

w = @(z) w_0*sqrt(1+(z/z_r).^2);
R_z = @(z) z./2./(z.^2+z_r^2);
phi_z = @(z) atan(z./z_r);
A = @(x,y,z) 1./w(z).*exp((-x.^2-y.^2)./w(z).^2);
Phi = @(x,y,z) k*((-x.^2-y.^2).*R_z(z)+z)-phi_z(z);

E_11 = @(x,y) A(x,y,0);
E_22 = @(x,y) A(x*cos(delta),y,x*sin(delta));
E_1 = @(x,y) E_11(x,y).*E_11(x,y);
E_2 = @(x,y) E_22(x,y).*E_22(x,y);
E_3 = @(x,y) E_11(x,y).*E_22(x,y).*exp(-i*(Phi(x,y,0)-Phi(x*cos(delta),y,x*sin(delta))));
Phi_2 = @(x,y) Phi(x.*cos(delta),y,x.*sin(delta));

x = -160:160;
y = -128:128;
mu = 30e-6;
%%
for j = 1:length(x)-1
    for k = 1:length(y)-1
        xmin = x(j)*mu;
        xmax = (x(j)+1)*mu;
        ymin = y(k)*mu;
        ymax = (y(k)+1)*mu;
        I_1(j,k) = integral2(E_1,xmin,xmax,ymin,ymax,'Method','iterated');
        I_2(j,k) = integral2(E_2,xmin,xmax,ymin,ymax,'Method','iterated');
        I_3(j,k) = integral2(E_3,xmin,xmax,ymin,ymax,'Method','iterated');
%         EE(j,k) = I_3^2/(I_1*I_2);
    end
    j
end
%%
mesh(angle(I_3)')
xlabel('x轴像素')
ylabel('y轴像素')
title('CCD平面相位分布')
axis([0 320 0 256])