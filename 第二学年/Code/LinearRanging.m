clc
clear

lambda = 1064e-9; % 波长
w_0 = 0.3125e-3; % 光斑半径

k = 2*pi/lambda; 
z_r = pi*1*w_0^2/lambda;

x = 0;
alpha = 0:1e-4:1e-2;

z = 5e-2;
y = 0;

for j = 1:length(x)
    xm = x(j)*cos(alpha)-z*sin(alpha);
    ym = y*alpha;
    zm = z*cos(alpha)+ x(j)*sin(alpha);
    
    A1 = x(j)*z*cos(2*alpha)+(x(j)^2-z^2)/2*sin(2*alpha);
    B1 = z^2*cos(alpha).^2+x(j)^2*sin(alpha).^2+x(j)*z*sin(2*alpha)+z_r^2;
    PA1 = -2*x(j)*z*sin(2*alpha)+(x(j)^2-z^2)*cos(2*alpha);
    PB1 = -z*sin(2*alpha)+x(j)^2*sin(2*alpha)+2*x(j)*z*cos(2*alpha);
    
    P1 = (PA1.*B1-PB1.*A1)./B1.^2;
    
    A2 = (zm.^2-z_r^2).*(xm.^2+ym.^2).*sin(alpha);
    B2 = zm.^3+zm*z_r;
    PA2 =  2*zm.*xm.*(zm.^2+ym.^2).*sin(alpha) + (-2*xm.*zm.*sin(alpha)+cos(alpha).*(xm.^2+ym.^2)).*(zm.^2-z_r^2);
    PB2 = 3*zm.^2.*xm+z_r*xm;
    
    P2 = 0.5*(PA2.*B2-PB2.*A2)./B2.^2;
    
    P3 = (z_r*cos(alpha).*(z_r^2+zm.^2)-z_r*sin(alpha).*(2.*(zm.*xm)))./(z_r^2+zm.^2).^2;
    
    a1 = A1./B1;
    a2 = A2./B2;
    a3 = sin(alpha)./z_r./(1+(zm/z_r).^2);
    
    a_alpha(:,j) = -k*(P1+P2)-k*cos(alpha);
    
end







































% clear
% clc
% syms x y z;
% syms alpha;
% syms lr(x,y,z,alpha);
% 
% xm = x*cos(alpha) - z*sin(alpha);
% ym = y;
% zm = z*cos(alpha) + x*sin(alpha);
% 
% lambda = 1064e-9; % 波长
% w_0 = 0.3125e-3; % 光斑半径
% 
% k = 2*pi/lambda; 
% z_r = pi*1*w_0^2/lambda;
% 
% Rz1 =  z./2./(z.^2+z_r^2);%等相面曲率半径
% phiz1 = atan(z./z_r);%高斯光的相位因子
% 
% Rz2 =  zm./2./(zm.^2+z_r^2);%等相面曲率半径
% phiz2 = atan(zm./z_r);%高斯光的相位因子
% 
% phi1 = k*((x.^2+y.^2)./2./Rz1+z_r) - phiz1;
% phi2 = k*((xm.^2+ym.^2)./.2./Rz2+z_r) - phiz2;
% 
% Phi = phi1-phi2;
% 
% a = diff(Phi,x);
% 
% lr = diff(a,alpha);
