clear
x = -10:0.1:10;
y = x;
theta = pi/3;
[X,Y] = meshgrid(x,y);
hold on
Z = -1/tan(theta)*X;
mesh(X,Y,Z)
%%
hold on
%法线
y1 = 0*x+1;
z1 = tan(theta)*x;
plot3(x,y1,z1)
hold on
grid on
xlim([-10,10])
ylim([-10,10])
zlim([-10,10])
axis equal
%%
AB = [-1,1,0];
n1 = [cos(theta),0,sin(theta)];
n2 = cross(AB,n1);
n3 = cross(n1,n2);
BC = -sum(AB.*n1)*n1/norm(n1)+sum(AB.*n3)*n3/norm(n3)^2;
z2 = BC(3)/BC(1)*x;
y2 = BC(2)/BC(1)*x+1;
%%
hold on 
plot3(x,y2,z2);
%%
y3 = -x+1;
z3 = 0*x;
hold on
plot3(x,y3,z3)
%% 
BC1 = [cos(2*theta),-1,sin(2*theta)];
y4 = 1/cos(2*theta)*x+1;
z4 = tan(2*theta)*x;
hold on
plot3(x,y4,z4)
legend('反射面','法线','反射光','入射光','反射光2')
