clear
clc
close
%��˹��Ļ�����������߰뾶�Ͳ���
w_0 = 0.3125e-3;
lambda = 1064e-9;
% ��������
k = 2*pi/lambda;
z_r = pi*1*w_0^2/lambda;

w = @(z) w_0*sqrt(1+(z/z_r).^2);%��˹��Ĺ�߰뾶
R_z = @(z) z./2./(z.^2+z_r^2);%���������ʰ뾶
phi_z = @(z) atan(z./z_r);%��˹�����λ����
A = @(x,y,z) 1./w(z).*exp((-x.^2-y.^2)./w(z).^2);%
Phi = @(x,y,z) k*((x.^2+y.^2).*R_z(z)+z)-phi_z(z);
%% ������
z = -1:1e-2:1;
theta = sqrt(lambda/pi/z_r);% Զ�������
z1 = theta*z;

x = linspace(-0.640e-3,0.640e-3,640);
wx = real(A(x,0*x,0*x).*exp(1i*Phi(x,0*x,0*x)));%��ǿ��x��ķֲ�

% ��λ�ֲ�

y = linspace(-0.512e-3,0.512e-3,512);
[X, Y] = meshgrid(x, y);
[n,h] = size(X);
z0 = 0;%10e-2;
phi = Phi(X, Y, z0*ones(n,h));%z = 5cm������λ�ֲ�
%%
% ������������
alpha = 1e-3;% ���ƫת�Ƕ�
beta = 0e-3;


E_L = A(X,Y,5e-2*ones(n,h)).*exp(1i*Phi(X,Y,z0*ones(n,h)));
Ralpha = [cos(alpha),0,sin(alpha);0,1,0;-sin(alpha),0,cos(alpha)];%��y��ת��
Rbeta = [1,0,0;0,cos(beta),sin(beta);0,-sin(beta),cos(beta)];%��x��ת��
for m = 1:n
    for j = 1:h
        R = [X(m,j) Y(m,j) z0]*Ralpha*Rbeta;
        X1(m,j) = R(1);Y1(m,j) = R(2);Z1(m,j) = R(3);
    end
end
phi_2 = Phi(X1,Y1,Z1);
E_R = A(X1,Y1,Z1).*exp(1i*Phi(X1,Y1,Z1));
E_Phi = Phi(X,Y,z0*ones(n,h)) - Phi(X1,Y1,Z1);

%����������С���˷����õ�����ϵ��
p = ParameterInMatrix(X,Y,E_Phi);
sprintf('����ƫת�Ƕ�Ϊ��alpha = %d, beta = %d;����ϵ���ֱ�Ϊ��a = %d, b = %d, c = %d',alpha, beta, p(1),p(2),p(3))
sprintf('�����󣬹���ƫת�Ƕ�Ϊ��alpha = %d, beta = %d;����õ���ƫת�Ƕ�Ϊ��alpha = %d, beta = %d',alpha, beta, p(1)/(-6074517), p(2)/(-6074517))
%% ��ͼ��
close all
subplot 121
plot(z,w(z),'k',z,-w(z),'k',z,z1,'k--',z,-z1,'k--',z,0*z,'k--')%��߰뾶��z��ı仯
xticks([0]);yticks([])
title('�����뾶��z��ı仯')

subplot 122
plot(x,wx,'k');%��ǿ��x��ķֲ�
xticks([-1e-2 0 1e-2]);yticks([])
title('��ǿ��x���ϵı仯')

figure
subplot 121
mesh(X*1e3,Y*1e3,phi)%��λ��̽�����ϵķֲ�����λ�����ֵ-2.949e5
title('������δƫתʱ����λ�ֲ��������䣩')
xlim([-0.64,0.64])
ylim([-0.512,0.512])
xlabel('x��/mm')
ylabel('y��/mm')
subplot 122
mesh(X*1e3,Y*1e3,phi_2)% ��б��������λ�ֲ�-2.935e5
title('������ƫתʱ����λ�ֲ��������䣩')
xlabel('x��/mm')
ylabel('y��/mm')
xlim([-0.64,0.64])
ylim([-0.512,0.512])

figure
subplot 121
mesh(X*1e3,Y*1e3,phi);
title('����δƫתʱ����λ�ֲ�(������)')
xlabel('x��/mm')
ylabel('y��/mm')
xlim([-0.64,0.64])
ylim([-0.512,0.512])

subplot 122
mesh(X*1e3,Y*1e3,angle(E_R));
title('����ƫתʱ����λ�ֲ��������䣩')
xlabel('x��/mm')
ylabel('y��/mm')
xlim([-0.64,0.64])
ylim([-0.512,0.512])

figure
subplot 121
mesh(X*1e3,Y*1e3,E_Phi)
title('�����źŵ���λ�ֲ�')
xlabel('x��/mm');
ylabel('y��/mm')
xlim([-0.64,0.64])
ylim([-0.512,0.512])

subplot 122
mesh(X*1e3,Y*1e3,angle(exp(1i*E_Phi)));
title('�����źŵ���λ�ֲ�������λ���䣩')
xlabel('x��/mm');
ylabel('y��/mm')
xlim([-0.64,0.64])
ylim([-0.512,0.512])
