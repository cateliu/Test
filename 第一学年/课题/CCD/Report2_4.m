%% CCD�ļ���
%����ĽǶ�ȷ������y��תdelta��alpha = 0��delta = 0.2mrad��
%CCD��΢Ԫ���Ϊ30e-6 X 30e-6������Ϊ 320x256��
% ������˹�����䣬���������ͬ�����ǶȲ�ͬ��
% ����Ƶ�ʲ��Ӱ��
clear
clc
format long
% delta = 0.1e-3;
f = 1623;% Ƶ�ʲ�
T = 1/f;

w_0 = 0.3125e-3;
lambda = 1064e-9;
k = 2*pi/lambda;
z_r = pi*1*w_0^2/lambda;

w = @(z) w_0*sqrt(1+(z/z_r).^2);
R_z = @(z) z./2./(z.^2+z_r^2);
phi_z = @(z) atan(z./z_r);
A = @(x,y,z) 1./w(z).*exp((-x.^2-y.^2)./w(z).^2);
Phi = @(x,y,z) k*((-x.^2-y.^2).*R_z(z)+z)-phi_z(z);

E_11 = @(x,y) A(x,y,0);
E_1 = @(x,y) E_11(x,y).*E_11(x,y);


Phi_2 = @(x,y) Phi(x.*cos(delta),y,x.*sin(delta));

x = -160:159;
y = -128:127;
mu = 30e-6;
theta = 1e-3:1e-3:50e-3;
%%  ����60������
for t = 1:length(theta)
%     t1 = T/40*t;
    delta = theta(t);
    E_22 = @(x,y) A(x*cos(delta),y,x*sin(delta));
    E_2 = @(x,y) E_22(x,y).*E_22(x,y);
    
    E_3 = @(x,y) E_11(x,y).*E_22(x,y).*exp(i*(Phi(x,y,0)-Phi(x*cos(delta),y,x*sin(delta))));
    I(:,:,t) = c(x,y,E_1,E_2,E_3,mu);
    t/length(theta)
end
%% ����ֱ����
temp = deletAverage(I);
%% ������λ
phi_init = atan(imag(temp(:,:,1))./real(temp(:,:,1)));
subplot 221
mesh(phi_init)
xlabel('x������');ylabel('y������');zlabel('��λ [rad]');title('�������λ')
%% ������λ
format long
phi_rev = FourPoint(I,10);
subplot 222
mesh(X,Y,phi_rev)
xlabel('x������');ylabel('y������');zlabel('��λ [rad]');title('ʹ���ĵ��㷨����õ�����λ�ֲ�');axis([-128 128 -160 160])
%%  ������λ����
phi_smoth = Smooth(phi_rev);
subplot 223
mesh(phi_smoth)
xlabel('x������');ylabel('y������');zlabel('��λ [rad]');title('������������λ�ֲ�')
[r,v] = find(temp(:,:,1) == max(max(temp(:,:,1)))); 
phi_smoth_fit = phi_smoth - (phi_smoth(r,v)-phi_rev(r,v));
subplot 224
mesh(phi_rev)
hold on
mesh(phi_smoth_fit)
xlabel('x������');ylabel('y������');zlabel('��λ [rad]');title('������������λ�ֲ���δ������ʼ��λ�ֲ��ıȽ�')
%%
function T1 = c(x,y,E_1,E_2,E_3,mu)
    for j = 1:length(x)
        for k = 1:length(y)
            xmin = x(j)*mu;
            xmax = (x(j)+1)*mu;
            ymin = y(k)*mu;
            ymax = (y(k)+1)*mu;
            I_1 = integral2(E_1,xmin,xmax,ymin,ymax,'Method','iterated');
            I_2 = integral2(E_2,xmin,xmax,ymin,ymax,'Method','iterated');
            I_3 = integral2(E_3,xmin,xmax,ymin,ymax,'Method','iterated');
            T1(j,k) = 2*I_3./sqrt(I_1.*I_2);
        end
        j
    end
end

function T = deletAverage(I)
% ��������I��ֱ������
[jmax, mmax, imax] = size(I);
for j = 1:jmax
    for m = 1:mmax
        T(j,m,:) = I(j,m,:)-mean(real(I(j,m,:)));
    end
end
end

function T = atanphi(I)
[jmax, mmax, imax] = size(I);
for j = 1:imax
    T(:,:,j) = atan(imag(I(:,:,j))./real(I(:,:,j)));
end
end

function T = FourPoint(I,deltak)
% �ĵ��㷨 �����I��һ��������λ
[jmax, mmax, imax] = size(I);
Re = real(I);
for j = 1:jmax
    for k = 1:mmax
    T(j,k)=atan((Re(j,k,1+3*deltak)-Re(j,k,1+deltak))/(Re(j,k,1)-Re(j,k,1+2*deltak)));
    end
end
end

function T = Smooth(I)
% ������λ����
[r,v] = size(I);
A_1 = I;
for j = 1:r
     T(j,:)=compent(A_1(j,:));
end
T(:,1)=compent(T(:,1));
for j = 1:r
     T(j,:)=compent(T(j,:));
end
end

function T = compent(C)
% ����һ�������������ݴ������䣬������������Ե�һ��Ԫ��Ϊ��׼������������ķ����õ�������ֵ
l = length(C);
delta = 0;%������λ��
delta1 = 0; %����ƽ����
T(1) = C(1);
for j = 1:l-1
    chazhi = C(j+1)-C(j);
        if abs(chazhi)>=3
        if j == 1
            delta1 = C(4)-C(3);
%         elseif j==2
%             delta1 = C(4)-C(3);
        else
            delta1 = C(j)-C(j-1);
        end
        delta = delta+chazhi-delta1;
        end
    T(j+1) = C(j+1)-delta;
end
end
