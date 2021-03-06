clear
clc

w_0 = 0.3125e-3;
lambda = 1064e-9;

k = 2*pi/lambda;
z_r = pi*1*w_0^2/lambda;

w = @(z) w_0*sqrt(1+(z/z_r).^2);
R_z = @(z) z./2./(z.^2+z_r^2);
phi_z = @(z) atan(z./z_r);
A = @(x,y,z) 1./w(z).*exp((-x.^2-y.^2)./w(z).^2);
Phi = @(x,y,z) k*((-x.^2-y.^2).*R_z(z)+z)-phi_z(z);

theta = -10e-3:1e-4:10e-3;
r = [5e-4,5e-5,5e-6];

E_1 = @(x,y) A(x,y,0).*A(x,y,0);
% O = zeros(length(r),length(theta));
for j = 1:length(theta)
    E_2 = @(x,y) A(x*cos(theta(j)),y,x*sin(theta(j))).*A(x*cos(theta(j)),y,x*sin(theta(j)));
    E_3 = @(x,y) A(x,y,0).*A(x*cos(theta(j)),y,x*sin(theta(j))).*exp(-1i*(Phi(x,y,0)-Phi(x*cos(theta(j)),y,x*sin(theta(j)))));
    for m = 1:length(r)
        ymax = @(x) sqrt(r(m)^2-x.^2);
        I1(j,m) = integral2(E_1,0,r(m),0,ymax);
        I2(j,m) = integral2(E_2,0,r(m),0,ymax);
        I3(j,m) = integral2(E_3,0,r(m),0,ymax);
        O(j,m) = abs(I3(j,m)/(I1(j,m)*I2(j,m))^0.5);
    end
%     sprintf(num2str(j/length(theta))+' \n')
end
Eta = O.^2;
% r = 1e6*r;
% plot(r,Eta(1,:),r,Eta(2,:),r,Eta(3,:))
% plot(theta*1e3,O)
