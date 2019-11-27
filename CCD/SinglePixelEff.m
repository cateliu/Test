clc
clear
format long
lambda = 1064e-9;
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
E = @(x,y,z) A(x,y,z).*exp(-i*Phi(x,y,z));
E_1 = @(x,y) A(x,y,0).^2;
I_1 = integral2(E_1,-mu,mu,-mu,mu,'Method','iterated');
% I_1 = Inter(E_1,500,500);
% I_1 = sum(sum(E_1(X,Y)))*ds^2;

alpha = -8e-3:1e-4*2:8e-3;
delta = alpha;


%%
for j = 1:length(alpha)
	for m = 1:length(delta)
%         X = @(x,y) [x*cos(alpha(j)), y*cos(delta(m)) - x*sin(alpha(j))*sin(delta(m)), y*sin(delta(m)) + x*cos(delta(m))*sin(alpha(j))];
		E_22 = @(x,y) A(x*cos(alpha(j)), y*cos(delta(m)) - x*sin(alpha(j))*sin(delta(m)), y*sin(delta(m)) + x*cos(delta(m))*sin(alpha(j)));
        E_2 = @(x,y) E_22(x,y).*E_22(x,y);
		E_3 = @(x,y) A(x,y,0).*E_22(x,y).*exp(-i*(Phi(x,y,0)-Phi(x*cos(alpha(j)), y*cos(delta(m)) - x*sin(alpha(j))*sin(delta(m)), y*sin(delta(m)) + x*cos(delta(m))*sin(alpha(j)))));
		I_2 = integral2(E_2,-mu,mu,-mu,mu,'Method','iterated');
		I_3 = integral2(E_3,-mu,mu,-mu,mu,'Method','iterated');
%         I_2 = sum(sum(E_2(X,Y)))*ds^2;
%         I_3 = sum(sum(E_3(X,Y)))*ds^2;
		Eta(j,m) = I_3^2/I_1/I_2;
    end
    j
end
[A,B]=meshgrid(alpha*1e3,delta*1e3);
mesh(A,B,real(Eta));
xlabel('Horizon tilt /mrad')
ylabel('Vertical /mrad')
zlabel('Hetreodyne efficiency')