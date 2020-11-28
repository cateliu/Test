clear
clc
close
format long
%��˹��Ļ�����������߰뾶�Ͳ���

w_0 = 1e-3 ;%+ 0.3125e-3*(d-1)*0.5;
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

% x = linspace(-3.20e-3,3.20e-3,320);
x = 60e-6*[[1:320]-160];
wx = real(A(x,0*x,0*x).*exp(1i*Phi(x,0*x,0*x)));%��ǿ��x��ķֲ�

% ��λ�ֲ�
y = 60e-6*[[1:256]-128];
% y = linspace(-2.56e-3,2.56e-3,256);
[X, Y] = meshgrid(x, y);
[n,h] = size(X);
%%
for f = 1:4
    theta1 = 1e-8;-1e-6:1e-7:1e-6;

    z0 = 0;10e-2*(f-1);%10e-2;
    phi = Phi(X, Y, z0*ones(n,h));%z = 5cm������λ�ֲ�

    for q = 1:length(theta1)
    % ������������
        alpha = theta1(q);% ���ƫת�Ƕ�
        beta = 0;10e-3;


        E_L = A(X,Y,z0*ones(n,h)).*exp(1i*Phi(X,Y,z0*ones(n,h)));
        Ralpha = [cos(alpha),0,sin(alpha);0,1,0;-sin(alpha),0,cos(alpha)];%��y��ת��
        Rbeta = [1,0,0;0,cos(beta),sin(beta);0,-sin(beta),cos(beta)];%��x��ת��
        for m = 1:n
            for j = 1:h
                R = [X(m,j) Y(m,j) z0]*Ralpha*Rbeta;
                X1(m,j) = R(1);Y1(m,j) = R(2);Z1(m,j) = R(3);
            end
        end
        phi_2 = Phi(X1,Y1,Z1);
        E_R = A(X1,Y1,Z1).*exp(1i*(Phi(X1,Y1,Z1)+pi/2*(f-1)));
        P(:,:,f) = (E_L+E_R).*conj(E_L+E_R);
        E_Phi = Phi(X,Y,z0*ones(n,h)) -( Phi(X1,Y1,Z1));
        phi(:,:,f) = E_Phi;
        %����������С���˷����õ�����ϵ��
        phase = PhaseUnwrapping(atan2(sin(E_Phi),cos(E_Phi)),3 , 2*pi);
        p = ParameterInMatrix(phase, 60e-6);%+ (d-1)*rand(512,640)*10^(-(f-1)*0.2)
        a(q, f) = p(1);
    end

end
%%
% a1 = a/(-k);
% semilogx(10.^(-((1:31)-1)*0.3),(a1-10e-3)*1e6,'LineWidth',2)
% xlim([10^(-5),1])
% xlabel("��λ���/rad")
% ylabel("������/urad")
% xticks([10^(-5) 10^(-3) 1])
% xticks([10^(-5) 10^(-3) 10^(-2) 1])
% legend("û����λ���","������λ���")
% grid on
figure
phi1 = atan((P(:,:,4)-P(:,:,2))./((P(:,:,1)-P(:,:,3))));
mesh(phi1*1e6);xlim([1 320]);ylim([1 256]);xlabel("x��");ylabel("y��")
figure
subplot 221
mesh(P(:,:,1))
xlim([0 320]);xticks([0 320])
ylim([0 256]);yticks([256])
title("t=0")
subplot 222
mesh(P(:,:,2))
xlim([0 320]);xticks([0 320])
ylim([0 256]);yticks([256])
title("t=T/4")
subplot 223
mesh(P(:,:,3))
xlim([0 320]);xticks([0 320])
ylim([0 256]);yticks([256])
title("t=T/2")
subplot 224
mesh(P(:,:,4));
xlim([0 320]);xticks([0 320])
ylim([0 256]);yticks([256])
title("t=3T/4")
%%
figure
[W,Q] = meshgrid(151:483,103:437);
mesh(W,Q, phi1(103:437,151:483));xlim([151 483]);ylim([103 437]);xlabel("x��");ylabel("y��")

figure
u = phi1(103:437,151:483);
plot(103:437,u(:,1))
xlim([103,437])
yticks([-0.02,0 0.02])
grid on
xlabel("y��")
ylabel("��λ/rad")
%%

Phi1 = PhaseUnwrapping(phi1,0.5,pi);
% for j = 1:512
%     Phi1(j,:) = PhaseUnwrapping(phi1(j,:),0.5,pi);
% end
% mesh(Phi1)
% plot((0.3125e-3 + 0.3125e-3*([1:20]-1)*0.2)*1e3, (a1-10e-3)*1e6,'LineWidth',2)
% xlim([0.3125 1.5])
% grid on
% yticks([-0.5 -0.15 0 0.5 1 1.5 2])
% xlabel("��߰뾶/mm")
% ylabel("������/urad")

% plot(theta1'*1e3,(a)*1e3,'LineWidth',2)
% plot(0.3125e-3 + 0.3125e-3*([1:20]-1)*0.2, (a1-10e-3)*1e3,'LineWidth',2)
% plot(theta1*1e3,a1)
% xticks([-100 -50 0 50 100])
% yticks([-100 -50 0 50 100])
% xlim([-100,100])
% xlabel("ƫת�Ƕ�/mrad")
% ylabel('����ֵ/mrad')
% grid on
% legend("w=w0","w=2w0","w=3w0","w=4w0","w=5w0","w=6w0","w=7w0","w=8w0")

%legend("z=0cm","z=20cm","z=40cm","z=60cm","z=80cm")
% axis equal
% sprintf('����ƫת�Ƕ�Ϊ��alpha = %d, beta = %d;����ϵ���ֱ�Ϊ��a = %d, b = %d, c = %d',alpha, beta, p(1),p(2),p(3))
% sprintf('�����󣬹���ƫת�Ƕ�Ϊ��alpha = %d, beta = %d;����õ���ƫת�Ƕ�Ϊ��alpha = %d, beta = %d',alpha, beta, p(1)/(-6074517), p(2)/(-6074517))
%% ��ͼ��
% close all
% subplot 121
% plot(z,w(z),'k',z,-w(z),'k',z,z1,'k--',z,-z1,'k--',z,0*z,'k--')%��߰뾶��z��ı仯
% xticks([0]);yticks([])
% title('�����뾶��z��ı仯')
% 
% subplot 122
% plot(x,wx,'k');%��ǿ��x��ķֲ�
% xticks([-1e-2 0 1e-2]);yticks([])
% title('��ǿ��x���ϵı仯')
% 
% figure
% subplot 121
% mesh(X*1e3,Y*1e3,phi)%��λ��̽�����ϵķֲ�����λ�����ֵ-2.949e5
% title('������δƫתʱ����λ�ֲ�')
% xlim([-6.4,6.4])
% ylim([-5.12,5.12])
% xlabel('x��/mm')
% ylabel('y��/mm')
% subplot 122
% mesh(X*1e3,Y*1e3,phi_2)% ��б��������λ�ֲ�-2.935e5
% title('������ƫתʱ����λ�ֲ�')
% xlabel('x��/mm')
% ylabel('y��/mm')
% xlim([-6.4,6.4])
% ylim([-5.12,5.12])
% 
% figure
% subplot 121
% mesh(X*1e3,Y*1e3,phi);
% title('����δƫתʱ����λ�ֲ�(������)')
% xlabel('x��/mm')
% ylabel('y��/mm')
% xlim([-6.4,6.4])
% ylim([-5.12,5.12])
% 
% subplot 122
% mesh(X*1e3,Y*1e3,angle(E_R));
% title('����ƫתʱ����λ�ֲ��������䣩')
% xlabel('x��/mm')
% ylabel('y��/mm')
% xlim([-6.4,6.4])
% ylim([-5.12,5.12])
% 
% figure
% subplot 121
% mesh(X*1e3,Y*1e3,E_Phi)
% title('�����źŵ���λ�ֲ�')
% xlabel('x��/mm');
% ylabel('y��/mm')
% xlim([-6.4,6.4])
% ylim([-5.12,5.12])
% 
% subplot 122
% mesh(X*1e3,Y*1e3,angle(exp(1i*E_Phi)));
% title('�����źŵ���λ�ֲ�������λ���䣩')
% xlabel('x��/mm');
% ylabel('y��/mm')
% xlim([-6.4,6.4])
% ylim([-5.12,5.12])
