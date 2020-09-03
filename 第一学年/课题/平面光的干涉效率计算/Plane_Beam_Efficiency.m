clear
clc
alpha = -1e-4:1e-6:1e-4;
delta = -1e-4:1e-6:1e-4;
lambda = 1064e-9;
k_1 = 2*pi/lambda;
k1 = k_1*[0 0 1];
%%
%平面光的干涉效率计算，干涉角在(-0.1mrad,1mrad)内选取，正方形边长分别为0.2mm,2mm,20mm，分别对于q,q1,q2
%%
%先绕x轴旋转alpha角，再绕y轴旋转delta角
for j = 1:length(alpha)
	for m = 1:length(delta)
        k2 = k1*rx(alpha(j))*ry(delta(m));
        f_xyz =@(x,y)  exp(i.*((k1(1)-k2(1)).*x+(k1(2)-k2(2)).*y));
        q(j,m)=(real(integral2(f_xyz,-1e-4,1e-4,-1e-4,1e-4,'Method','iterated')))^2/(2e-4)^4;% 外差分布
        q1(j,m) = (real(integral2(f_xyz,-0.001,0.001,-0.001,0.001,'Method','iterated')))^2/(4e-6^2);
        q2(j,m) = (real(integral2(f_xyz,-0.01,0.01,-0.01,0.01,'Method','iterated')))^2/(4e-4^2);
        
    end
    fprintf('\n %d \n',j)
end

%%
[X,Y]=meshgrid(delta*1e3,alpha*1e3);
figure
mesh(X,Y,q)
xlabel('horizontal tilt /mrad','FontName','Times New Roman','FontSize',12);
 ylabel('vertical tilt /mrad','FontName','Times New Roman','FontSize',12);
 zlabel('heterodyne efficiency','FontName','Times New Roman','FontSize',12);
 title("外差效率随干涉角度的变化-0.2mm",'FontName','宋体','FontSize',12);

figure
 mesh(X,Y,q1)
xlabel('horizontal tilt /mrad','FontName','Times New Roman','FontSize',12);
 ylabel('vertical tilt /mrad','FontName','Times New Roman','FontSize',12);
  zlabel('heterodyne efficiency','FontName','Times New Roman','FontSize',12);
 title("外差效率随干涉角度的变化-2mm",'FontName','宋体','FontSize',12);
%   
figure
 mesh(X,Y,q2)
xlabel('horizontal tilt /mrad','FontName','Times New Roman','FontSize',12);
 ylabel('vertical tilt /mrad','FontName','Times New Roman','FontSize',12);
  zlabel('heterodyne efficiency','FontName','Times New Roman','FontSize',12);
 title("外差效率随干涉角度的变化-20mm",'FontName','宋体','FontSize',12);
 
