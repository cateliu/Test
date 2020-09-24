clear
clc
close
%高斯光的基本参数，光斑半径和波长
w_0 = 0.3125e-3;
lambda = 1064e-9;
% 衍生参数
k = 2*pi/lambda;
z_r = pi*1*w_0^2/lambda;

w = @(z) w_0*sqrt(1+(z/z_r).^2);%高斯光的光斑半径
R_z = @(z) z./2./(z.^2+z_r^2);%等相面曲率半径
phi_z = @(z) atan(z./z_r);%高斯光的相位因子
A = @(x,y,z) 1./w(z).*exp((-x.^2-y.^2)./w(z).^2);%
Phi = @(x,y,z) k*((x.^2+y.^2).*R_z(z)+z)-phi_z(z);
%% 数据区
z = -1:1e-2:1;
theta = sqrt(lambda/pi/z_r);% 远场发射角
z1 = theta*z;

x = linspace(-0.640e-3,0.640e-3,640);
wx = real(A(x,0*x,0*x).*exp(1i*Phi(x,0*x,0*x)));%光强随x轴的分布

% 相位分布

y = linspace(-0.512e-3,0.512e-3,512);
[X, Y] = meshgrid(x, y);
[n,h] = size(X);
z0 = 0;%10e-2;
phi = Phi(X, Y, z0*ones(n,h));%z = 5cm处的相位分布
%%
% 光束发生干涉
alpha = 1e-3;% 最大偏转角度
beta = 0e-3;


E_L = A(X,Y,5e-2*ones(n,h)).*exp(1i*Phi(X,Y,z0*ones(n,h)));
Ralpha = [cos(alpha),0,sin(alpha);0,1,0;-sin(alpha),0,cos(alpha)];%绕y轴转动
Rbeta = [1,0,0;0,cos(beta),sin(beta);0,-sin(beta),cos(beta)];%绕x轴转动
for m = 1:n
    for j = 1:h
        R = [X(m,j) Y(m,j) z0]*Ralpha*Rbeta;
        X1(m,j) = R(1);Y1(m,j) = R(2);Z1(m,j) = R(3);
    end
end
phi_2 = Phi(X1,Y1,Z1);
E_R = A(X1,Y1,Z1).*exp(1i*Phi(X1,Y1,Z1));
E_Phi = Phi(X,Y,z0*ones(n,h)) - Phi(X1,Y1,Z1);

%对曲面做最小二乘法，得到三个系数
p = ParameterInMatrix(X,Y,E_Phi);
sprintf('光束偏转角度为：alpha = %d, beta = %d;三个系数分别为：a = %d, b = %d, c = %d',alpha, beta, p(1),p(2),p(3))
sprintf('修正后，光束偏转角度为：alpha = %d, beta = %d;计算得到的偏转角度为：alpha = %d, beta = %d',alpha, beta, p(1)/(-6074517), p(2)/(-6074517))
%% 绘图区
close all
subplot 121
plot(z,w(z),'k',z,-w(z),'k',z,z1,'k--',z,-z1,'k--',z,0*z,'k--')%光斑半径随z轴的变化
xticks([0]);yticks([])
title('光束半径在z轴的变化')

subplot 122
plot(x,wx,'k');%光强随x轴的分布
xticks([-1e-2 0 1e-2]);yticks([])
title('光强在x轴上的变化')

figure
subplot 121
mesh(X*1e3,Y*1e3,phi)%相位在探测面上的分布，相位是相对值-2.949e5
title('光束在未偏转时的相位分布（无跳变）')
xlim([-0.64,0.64])
ylim([-0.512,0.512])
xlabel('x轴/mm')
ylabel('y轴/mm')
subplot 122
mesh(X*1e3,Y*1e3,phi_2)% 倾斜光束的相位分布-2.935e5
title('光束在偏转时的相位分布（无跳变）')
xlabel('x轴/mm')
ylabel('y轴/mm')
xlim([-0.64,0.64])
ylim([-0.512,0.512])

figure
subplot 121
mesh(X*1e3,Y*1e3,phi);
title('光束未偏转时的相位分布(有跳变)')
xlabel('x轴/mm')
ylabel('y轴/mm')
xlim([-0.64,0.64])
ylim([-0.512,0.512])

subplot 122
mesh(X*1e3,Y*1e3,angle(E_R));
title('光束偏转时的相位分布（有跳变）')
xlabel('x轴/mm')
ylabel('y轴/mm')
xlim([-0.64,0.64])
ylim([-0.512,0.512])

figure
subplot 121
mesh(X*1e3,Y*1e3,E_Phi)
title('干涉信号的相位分布')
xlabel('x轴/mm');
ylabel('y轴/mm')
xlim([-0.64,0.64])
ylim([-0.512,0.512])

subplot 122
mesh(X*1e3,Y*1e3,angle(exp(1i*E_Phi)));
title('干涉信号的相位分布（有相位跳变）')
xlabel('x轴/mm');
ylabel('y轴/mm')
xlim([-0.64,0.64])
ylim([-0.512,0.512])
