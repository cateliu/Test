% 这个脚本用来计算CCD探测面所有像元的探测效率在光的偏转大于1mrad时的干涉效率
% CCD的像元尺寸为30微米x30微米。阵列为256x320.
clear
clc
% CCD的尺寸
mu = 30e-6;
% 高斯光的参数
lambda = 1064e-9;
w_0 = 2e-3;
z_r = pi*w_0^2/lambda;
k = 2*pi/lambda;
% 基模高斯光的函数定义
w = @(z) w_0*sqrt(1+(z/z_r).^2);
phi_z = @(z) atan(z./z_r);
R_z = @(z) z./2./(z.^2+z_r^2);
A = @(x,y,z) 1./w(z).*exp((-x.^2-y.^2)./w(z).^2);
Phi = @(x,y,z) k*((-x.^2-y.^2).*R_z(z)+z)-phi_z(z);
E = @(x,y,z) A(x,y,z).*exp(-i*Phi(x,y,z));
% CCD的像素分布
y = -128:127;
x = -160:159;
%绕x轴旋转角度

for t = 1:10
    alpha = 1e-3*t;
    E_1 = @(x,y) A(x,y,0).*A(x,y,0);
    E_2 = @(x,y) A(x, y*cos(alpha), y*sin(alpha)).*A(x, y*cos(alpha), y*sin(alpha));
    E_3 = @(x,y) A(x,y,0).*A(x, y*cos(alpha), y*sin(alpha)).*exp(-i.*(Phi(x,y,0)-Phi(x, y*cos(alpha), y*sin(alpha))));
% 计算CCD的外差效率分布
    for j = 1:length(x)
        for l = 1:length(y)
            xmin = x(j)*mu;
            xmax = (x(j)+1)*mu;
            ymin = y(l)*mu;
            ymax = (y(l)+1)*mu;
            O_1(l,j) = integral2(E_1,xmin,xmax,ymin,ymax);
            O_2(l,j) = integral2(E_2,xmin,xmax,ymin,ymax);
            O_3(l,j) = integral2(E_3,xmin,xmax,ymin,ymax);                 
            Eta(l,j) = O_3(l,j)^2/O_1(l,j)/O_2(l,j);
        end
    end
    [X,Y] = meshgrid(x,y);
    mesh(X,Y,abs(Eta));
    xlabel('x轴像素')
    ylabel('y轴像素')
    axis([-160 160 -128 128 0 1])
    title("偏转角度为："+num2str(t)+"mrad")
    f = gcf;
    f.Color = 'w';
    frame = getframe(gcf);
    imind = frame2im(frame);
    [imind,cm] = rgb2ind(imind,256);
    if t==1
         imwrite(imind,cm,'CCD_Distr1.gif','gif', 'Loopcount',inf,'DelayTime',0.01);
    else
         imwrite(imind,cm,'CCD_Distr1.gif','gif','WriteMode','append','DelayTime',0.01);
    end
end

%%

% axis([-128 128 -160 160])