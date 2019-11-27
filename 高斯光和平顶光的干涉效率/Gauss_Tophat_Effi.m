%%  平顶光与高斯光的干涉效率计算
clc
clear
alpha =  -8e-3:1e-4:8e-3;
delta =  alpha;
w_0 = 0.3125e-3;

lambda = 1064e-9;
k = 2*pi/lambda;
z_r = pi*1*w_0^2/lambda;
w = @(z) w_0*sqrt(1+(z/z_r).^2);
R_z = @(z) z./2./(z.^2+z_r^2);
phi_z = @(z) atan(z./z_r);
E = @(x,y,z) 1./w(z).*exp((-x.^2-y.^2)./w(z).^2);
Phi = @(x,y,z) k*((-x.^2-y.^2).*R_z(z)+z)-phi_z(z);

r = 0.5e-3;
y = @(x) sqrt(0.5e-3^2-x.^2);
E_1_test = @(x,y) 0*x+1;
I_1 = integral2(E_1_test,-r,0,0,y,'Method','iterated');
%%
for j = 1:length(alpha)
    for m = 1:length(delta)
        E_2_test = @(x,y) exp(2*(-x.^2-y.^2)/w_0^2);
        E_3_test = @(x,y) exp((-x.^2-y.^2)/w_0^2).*exp(-i*k*(x*alpha(j)+y*delta(m)));
        
        I_2 = integral2(E_2_test,-r,0,0,y,'Method','iterated');
        I_3 = integral2(E_3_test,-r,0,0,y,'Method','iterated');
        O(j,m) = I_3^2/(I_1*I_2);
    end
    j
end
[x,y]=meshgrid(delta,alpha);
x = x*1e3;
y = y*1e3;
mesh(y,x,log10(abs(O))')
xlabel('horizontal tilt / mrad');
ylabel('vertical tilt / mrad');
hold on
contour(y,x,log10(abs(O))')