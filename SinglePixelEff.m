clc
clear
lambda = 1.064e-9;
w_0 = 0.3125e-3;
mu = 30e-6/2;
x  = -mu:1e-7:mu;
ds = x(2) - x(1);
[X,Y]=meshgrid(x);

k = 2*pi/lambda;
z_r = pi*1*w_0^2/lambda;

w = @(z) w_0*sqrt(1+(z/z_r).^2);
phi_z = @(z) atan(z./z_r);
R_z = @(z) z./2./(z.^2+z_r^2);
A = @(x,y,z) 1./w(z).*exp((-x.^2-y.^2)./w(z).^2);
Phi = @(x,y,z) k*((-x.^2-y.^2).*R_z(z)+z)-phi_z(z);

A_1 = @(x,y,z) 1./w(z).^2.*exp(-2*(x.^2+y.^2)./w(z).^2);
E = @(x,y,z) A(x,y,z).*exp(-i*Phi(x,y,z));
E_1 = @(x,y) A_1(x,y,0);
I_1 = integral2(E_1,-mu,mu,-mu,mu,'Method','iterated');
% I_1 = Inter(E_1,500,500);
I_1 = sum(sum(E_1(X,Y)))*ds^2;

alpha = -1e-4:1e-6:1e-4;
delta = alpha;


%%
for j = 1:length(alpha)
	for m = 1:length(delta)
%         X = @(x,y) [x*cos(alpha(j)), y*cos(delta(m)) - x*sin(alpha(j))*sin(delta(m)), y*sin(delta(m)) + x*cos(delta(m))*sin(alpha(j))];
		E_2 = @(x,y) A_1(x*cos(alpha(j)), y*cos(delta(m)) - x*sin(alpha(j))*sin(delta(m)), y*sin(delta(m)) + x*cos(delta(m))*sin(alpha(j)));
% 		E_3 = @(x,y) A(x,y,0).*A(x*cos(alpha(j)), y*cos(delta(m)) - x*sin(alpha(j))*sin(delta(m)), y*sin(delta(m)) + x*cos(delta(m))*sin(alpha(j))).*exp(-i*(Phi(x,y,0)-Phi(x*cos(alpha(j)), y*cos(delta(m)) - x*sin(alpha(j))*sin(delta(m)), y*sin(delta(m)) + x*cos(delta(m))*sin(alpha(j)))));
        E_m = @(x,y) E(x*cos(alpha(j)), y*cos(delta(m)) - x*sin(alpha(j))*sin(delta(m)), y*sin(delta(m)) + x*cos(delta(m))*sin(alpha(j)))+E(x,y,0);
        P = @(x,y) E_m(x,y).*conj(E_m(x,y));
        E_3 = @(x,y) (P(x,y)-E_1(x,y)-E_2(x,y))/2;
%         I_2(j,m) = Inter(E_2,500,500);
%         I_3 = Inter(E_3,500,500);
% 		I_2(j,m) = integral2(E_2,-mu,mu,-mu,mu,'Method','iterated');
% 		I_3 = integral2(E_3,-mu,mu,-mu,mu,'Method','iterated');
        I_2 = sum(sum(E_2(X,Y)))*ds^2;
        I_3 = sum(sum(E_3(X,Y)))*ds^2;
		Eta(j,m) = I_3^2/I_1/I_2;
    end
    j
end
[A,B]=meshgrid(alpha,delta);
mesh(A,B,abs(Eta));