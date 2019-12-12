% ����ű���������CCD̽����������Ԫ��̽��Ч���ڹ��ƫת����1mradʱ�ĸ���Ч��
% CCD����Ԫ�ߴ�Ϊ30΢��x30΢�ס�����Ϊ256x320.
clear
clc
% CCD�ĳߴ�
mu = 30e-6;

% ��˹��Ĳ���
lambda = 1064e-9;
w_0 = 0.3125e-3;
z_r = 2*pi*w_0^2/lambda;
k = 2*pi/lambda;
% ��ģ��˹��ĺ�������
w = @(z) w_0*sqrt(1+(z/z_r).^2);
phi_z = @(z) atan(z./z_r);
R_z = @(z) z./2./(z.^2+z_r^2);
A = @(x,y,z) 1./w(z).*exp((-x.^2-y.^2)./w(z).^2);
Phi = @(x,y,z) k*((-x.^2-y.^2).*R_z(z)+z)-phi_z(z);
E = @(x,y,z) A(x,y,z).*exp(-i*Phi(x,y,z));
% CCD�����طֲ�
x = -128:127;
y = -160:159;
%��x����ת�Ƕ�
alpha = 1e-3;
E_1 = @(x,y) A(x,y,0).*A(x,y,0);
E_2 = @(x,y) A(x, y*cos(alpha), y*sin(alpha)).*A(x, y*cos(alpha), y*sin(alpha));
E_3 = @(x,y) A(x,y,0).*A(x, y*cos(alpha), y*sin(alpha)).*exp(-i.*(Phi(x,y,0)-Phi(x, y*cos(alpha), y*sin(alpha))));


% ����CCD�����Ч�ʷֲ�
for j = 1:length(x)
    for l = 1:length(y)
        xmin = x(j)*mu;
        xmax = (x(j)+1)*mu;
        ymin = y(j)*mu;
        ymax = (y(j)+1)*mu;
        O_1 = integral2(E_1,xmin,xmax,ymin,ymax);
        O_2 = integral2(E_2,xmin,xmax,ymin,ymax);
        O_3 = integral2(E_3,xmin,xmax,ymin,ymax);
        Eta(l,j) = O_3^2/O_1/O_2;
    end
end
%%
[X,Y] = meshgrid(x,y);
mesh(X,Y,abs(Eta))
% axis([-128 128 -160 160])