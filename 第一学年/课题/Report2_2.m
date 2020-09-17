% 年终总结2.2 高斯光与高斯光在矩形面积下的外差效率分布
% 高斯光的参数是 lambda = 1064nm, z = 0, 矩形长度在1mm,0.1mm,0.01mm,
% 偏转角度为-8mrad-8mrad，w_0 = 0.3125mm
clear
clc
% 入射角度
alpha1 =  -8e-3:1e-5*2:8e-3;
delta1 =  0
% 矩形长度
l = [0.5e-3,0.05e-3,0.005e-3];
% 基模高斯光的定义
% 基本参数
w_0 = 0.3125e-3;
lambda = 1064e-9;
% 衍生参数
k = 2*pi/lambda;
z_r = pi*1*w_0^2/lambda;

w = @(z) w_0*sqrt(1+(z/z_r).^2);
R_z = @(z) z./2./(z.^2+z_r^2);
phi_z = @(z) atan(z./z_r);
A = @(x,y,z) 1./w(z).*exp((-x.^2-y.^2)./w(z).^2);
Phi = @(x,y,z) k*((-x.^2-y.^2).*R_z(z)+z)-phi_z(z);

% 设置积分参数
E_1 = @(x,y) A(x,y,0).*A(x,y,0);
I_1 = integral2(E_1,-l(1),l(1),-l(1),l(1));
O = zeros(801,801);
parfor j = 1:length(alpha1)
    v = zeros(1,length(delta1));
    for m =1:length(delta1)
        alpha = alpha1(j);
        delta = delta1(m);
        A2 = @(x,y) A(x*cos(delta) - y*sin(alpha)*sin(delta), y*cos(alpha), x*sin(delta) + y*cos(delta)*sin(alpha));
        E_2 = @(x,y) A2(x,y).*A2(x,y);
        E_3 = @(x,y) A(x,y,0).*A2(x,y).*exp(-i*(Phi(x,y,0)-Phi(x*cos(delta) - y*sin(alpha)*sin(delta), y*cos(alpha), x*sin(delta) + y*cos(delta)*sin(alpha))));
        I_2 = integral2(E_2,-l(1),l(1),-l(1),l(1));
        I_3 = integral2(E_3,-l(1),l(1),-l(1),l(1));
        v(m)= I_3^2/I_1/I_2;
    end
    O(j,:) = v;
    fprintf(num2str(j/length(delta1)*100)+" finished \n")
end
[X,Y] = meshgrid(delta1*1000,alpha1*1000);
meshc(X,Y,abs(O))
xlabel('horizon tilt [mrad]')
ylabel('vertical tilt [mrad]')
zlabel('Heterodyne Efficiency')
axis([-8 8 -8 8 0 1])












