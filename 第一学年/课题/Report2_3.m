% 年终总结2.3 CCD像素和QPD象限在同一个偏转下的干涉效率
% CCD像元为30微米，QPD半径为0.5mm
clear
clc
% 入射角度
alpha = -5e-3:1e-4:5e-3;
% 尺寸
mu = 30e-6;
r = 0.5e-3;
% 高斯光
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
I_1 = QPD(E_1,r);
C_1 = CCD(E_1,mu);
parfor j = 1:length(alpha)
    delta = alpha(j);
    E_2 = @(x,y) A(x*cos(delta),y,x*sin(delta)).*A(x*cos(delta),y,x*sin(delta));
    E_3 = @(x,y) A(x,y,0).*A(x*cos(delta),y,x*sin(delta)).*exp(-i*(Phi(x,y,0)-Phi(x*cos(delta),y,x*sin(delta))));
    I_2 = QPD(E_2,r);
    I_3 = QPD(E_3,r);
    C_2 = CCD(E_2,mu);
    C_3 = CCD(E_3,mu);
    I(j,:) = I_3.^2./I_1./I_2;
    C(j,:) = C_3.^2./C_2./C_1;
end
alpha = alpha*1e3;
%%
plot(alpha,abs(I(:,1)),'--',alpha,abs(I(:,2))-0.005,'--',alpha,abs(I(:,3))-0.01,'--',alpha,abs(I(:,4))-0.015,'--');
hold on
plot(alpha,abs(C(:,1)),alpha,abs(C(:,2))-0.005,alpha,abs(C(:,3))-0.01,alpha,abs(C(:,4))-0.015);
legend('A象限','B象限','C象限','D象限','CCD(1)','CCD(2)','CCD(3)','CCD(4)')
xlabel('偏转角度 [mrad]')
ylabel('外差效率')
% plot(alpha,angle(EtaA),alpha,angle(EtaB),alpha,angle(EtaC),alpha,angle(EtaD))

function T= QPD(E_1,r)
ymax = @(x) sqrt(r^2-x.^2);
ymin = @(x) -sqrt(r^2-x.^2);
T1 = integral2(E_1,-r,0,0,ymax,'Method','iterated');
T2 = integral2(E_1,0,r,0,ymax,'Method','iterated');
T3 = integral2(E_1,-r,0,0,ymin,'Method','iterated');
T4 = integral2(E_1,0,r,0,ymin,'Method','iterated');
T = [T1,T2,T3,T4];
end

function T = CCD(E_1,mu)
T1 = integral2(E_1,-mu,0,0,mu,'Method','iterated');
T2 = integral2(E_1,0,mu,0,mu,'Method','iterated');
T3 = integral2(E_1,-mu,0,0,-mu,'Method','iterated');
T4 = integral2(E_1,0,mu,0,-mu,'Method','iterated');
T = [T1,T2,T3,T4];
end