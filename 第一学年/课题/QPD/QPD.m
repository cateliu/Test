clc
clear
% w0=2mm,5x5mm square, lambda = 633nm
format long
lambda = 633e-9;
w_0 = 2e-3;
width  = 2.5e-3;
h = 2.5e-3;
g = 50e-6;

k = 2*pi/lambda;
z_r = pi*1*w_0^2/lambda;

w = @(z) w_0*sqrt(1+(z/z_r).^2);
phi_z = @(z) atan(z./z_r);
R_z = @(z) z./2./(z.^2+z_r^2);
A = @(x,y,z) 1./w(z).*exp((-x.^2-y.^2)./w(z).^2);
Phi = @(x,y,z) k*((-x.^2-y.^2).*R_z(z)+z)-phi_z(z);
E = @(x,y,z) A(x,y,z).*exp(-i*Phi(x,y,z));
E_1 = @(x,y) A(x,y,0).^2;

alpha = -200e-6:1e-6:200e-6;
QPD_A = integral2(E_1,-width,-g,g,h,'Method','iterated');
QPD_B = integral2(E_1,g,width,g,h,'Method','iterated');
QPD_C = integral2(E_1,-width,-g,-h,-g,'Method','iterated');
QPD_D = integral2(E_1,g,width,-h,-g,'Method','iterated');
%%
for j = 1:length(alpha)
    E_22 = @(x,y) A(x*cos(alpha(j)), y, -x*sin(alpha(j)));
    E_2 = @(x,y) E_22(x,y).*E_22(x,y);
    E_3 = @(x,y) A(x,y,0).*E_22(x,y).*exp(-i*(Phi(x,y,0)-Phi(x*cos(alpha(j)), y, -x*sin(alpha(j)))));
    QPD_A_2 = integral2(E_2,-width,-g,g,h,'Method','iterated');
    QPD_B_2 = integral2(E_2,g,width,g,h,'Method','iterated');
    QPD_C_2 = integral2(E_2,-width,-g,-h,-g,'Method','iterated');
    QPD_D_2 = integral2(E_2,g,width,-h,-g,'Method','iterated');
    
    QPD_A_3 = integral2(E_3,-width,-g,g,h,'Method','iterated');
    QPD_B_3 = integral2(E_3,g,width,g,h,'Method','iterated');
    QPD_C_3 = integral2(E_3,-width,-g,-h,-g,'Method','iterated');
    QPD_D_3 = integral2(E_3,g,width,-h,-g,'Method','iterated');
    I_A(j) = (QPD_A+QPD_A_2+2*QPD_A_3)/(QPD_A+QPD_A_2);
    I_B(j) = (QPD_B+QPD_B_2+2*QPD_B_3)/(QPD_B+QPD_B_2);
    I_C(j) = (QPD_C+QPD_C_2+2*QPD_C_3)/(QPD_C+QPD_C_2);
    I_D(j) = (QPD_D+QPD_D_2+2*QPD_D_3)/(QPD_D+QPD_D_2);
    Eta_A(j) = QPD_A_3/sqrt(QPD_A_2*QPD_A);
    Eta_B(j) = QPD_B_3/sqrt(QPD_B_2*QPD_B);
    Eta_C(j) = QPD_C_3/sqrt(QPD_C_2*QPD_C);
    Eta_D(j) = QPD_D_3/sqrt(QPD_D_2*QPD_D);
end
Angle=angle([Eta_A;Eta_B;Eta_C;Eta_D]);

pitch = Angle(1,:)+Angle(2,:)-Angle(3,:)-Angle(4,:);
yaw = Angle(1,:)+Angle(3,:)-Angle(2,:)-Angle(4,:);
plot(alpha,yaw)
figure 
plot(alpha,Angle(1,:),'*',alpha,Angle(2,:),'--',alpha,Angle(3,:),alpha,Angle(4,:),'D')
legend('A','B','C','D')
figure
plot(alpha,abs(I_A),alpha,abs(I_B),alpha,abs(I_C),alpha,abs(I_D))
figure 
plot(alpha,abs(Eta_A).^2,alpha,abs(Eta_B).^2,alpha,abs(Eta_C).^2,alpha,abs(Eta_D).^2)