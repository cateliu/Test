clc
clear
format long
lambda = 663e-9;
w_0 = 2e-3;
mu = 2.5e-3;30e-6/2;

k = 2*pi/lambda;
z_r = pi*1*w_0^2/lambda;

w = @(z) w_0*sqrt(1+(z/z_r).^2);
phi_z = @(z) atan(z./z_r);
R_z = @(z) z./2./(z.^2+z_r^2);
A = @(x,y,z) 1./w(z).*exp((-x.^2-y.^2)./w(z).^2);
Phi = @(x,y,z) k*((-x.^2-y.^2).*R_z(z)+z)-phi_z(z);
E = @(x,y,z) A(x,y,z).*exp(-i*Phi(x,y,z));
E_1 = @(x,y) A(x,y,0).^2;
I_1 = integral2(E_1,-mu,mu,-mu,mu,'Method','iterated');

alpha = -200e-6:1e-6:200e-6;
delta = 0;

y = @(x) sqrt(0.5e-3-x.^2);
Q_1 = integral2(E_1,-0.5e-3,0,0,y)
%%
for j = 1:length(alpha)
    E_22 = @(x,y) A(x*cos(alpha(j)), y, x*sin(alpha(j)));
    E_2 = @(x,y) E_22(x,y).*E_22(x,y);
    E_3 = @(x,y) A(x,y,0).*E_22(x,y).*exp(-i*(Phi(x,y,0)-Phi(x*cos(alpha(j)), y, x*sin(alpha(j)))));
    I_2 = integral2(E_2,-mu,mu,-mu,mu,'Method','iterated');
    I_3 = integral2(E_3,-mu,mu,-mu,mu,'Method','iterated');
    Eta(j) = I_3^2/I_1/I_2;
    Q_2 = integral2(E_2,-0.5e-3,0,0,y,'Method','iterated');
    Q_3 = integral2(E_3,-0.5e-3,0,0,y,'Method','iterated');
    eta(j) = Q_3^2/Q_1/Q_2;
    j
end
%%
plot(alpha,abs(eta))
hold on
plot(alpha,abs(Eta))
axis([-1e-2 1e-2 0 1])
xlabel('沿y轴偏转角')
ylabel('外差效率')
legend('QPD','CCD')

%%  差分信号的提取
QPD_A = integral2(E_1,-mu,0,0,mu,'Method','iterated');
QPD_B = integral2(E_1,0,mu,0,mu,'Method','iterated');
QPD_C = integral2(E_1,-mu,0,0,-mu,'Method','iterated');
QPD_D = integral2(E_1,0,mu,0,-mu,'Method','iterated');
for j = 1:length(alpha)
    E_22 = @(x,y) A(x*cos(alpha(j)), y, x*sin(alpha(j)));
    E_2 = @(x,y) E_22(x,y).*E_22(x,y);
    E_3 = @(x,y) A(x,y,0).*E_22(x,y).*exp(-i*(Phi(x,y,0)-Phi(x*cos(alpha(j)), y, x*sin(alpha(j)))));
    QPD_A_2 = integral2(E_2,-mu,0,0,mu,'Method','iterated');
    QPD_B_2 = integral2(E_2,0,mu,0,mu,'Method','iterated');
    QPD_C_2 = integral2(E_2,-mu,0,0,-mu,'Method','iterated');
    QPD_D_2 = integral2(E_2,0,mu,0,-mu,'Method','iterated');
    
    QPD_A_3 = integral2(E_3,-mu,0,0,mu,'Method','iterated');
    QPD_B_3 = integral2(E_3,0,mu,0,mu,'Method','iterated');
    QPD_C_3 = integral2(E_3,-mu,0,0,-mu,'Method','iterated');
    QPD_D_3 = integral2(E_3,0,mu,0,-mu,'Method','iterated');
    
    Eta_A(j) = QPD_A_3^2/QPD_A_2/QPD_A;
    Eta_B(j) = QPD_B_3^2/QPD_B_2/QPD_B;
    Eta_C(j) = QPD_C_3^2/QPD_C_2/QPD_C;
    Eta_D(j) = QPD_D_3^2/QPD_D_2/QPD_D;
end
Angle_A = angle(Eta_A)/2;
Angle_B = angle(Eta_B)/2;
Angle_C = angle(Eta_C)/2;
Angle_D = angle(Eta_D)/2;
pitch = Angle_A+Angle_B-Angle_C-Angle_D;
yaw = Angle_A+Angle_C-Angle_B-Angle_D;
plot(alpha,yaw)
figure
plot(alpha,Angle_A,'*',alpha,Angle_B,'+',alpha,Angle_C,alpha,Angle_D)
legend('A','B','C','D')