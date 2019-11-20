%% 显示H1和L1
clear
clc
tevent = 1126259462.422;
deltat = 5;

strain_H = h5read('H-H1_LOSC_4_V1-1126259446-32.hdf5','/strain/Strain');
gpsStart_H = cast(h5read('H-H1_LOSC_4_V1-1126259446-32.hdf5','/meta/GPSstart'),'double');
qmask_H = h5read('H-H1_LOSC_4_V1-1126259446-32.hdf5','/quality/simple/DQmask');
ts_H = h5readatt('H-H1_LOSC_4_V1-1126259446-32.hdf5','/strain/Strain','Xspacing');

strain_L = h5read('L-L1_LOSC_4_V1-1126259446-32.hdf5','/strain/Strain');
gpsStart_L = cast(h5read('L-L1_LOSC_4_V1-1126259446-32.hdf5','/meta/GPSstart'),'double');
qmask_L = h5read('L-L1_LOSC_4_V1-1126259446-32.hdf5','/quality/simple/DQmask');
ts_L = h5readatt('L-L1_LOSC_4_V1-1126259446-32.hdf5','/strain/Strain','Xspacing');

Time_H = gpsStart_H:ts_H:gpsStart_H+length(qmask_H);
index_H = find(Time_H>=tevent-deltat & Time_H<=tevent+deltat);
Time_L = gpsStart_L:ts_L:gpsStart_L+length(qmask_L);
index_L = find(Time_L>=tevent-deltat & Time_L<=tevent+deltat);

plot(Time_H(index_H)-tevent-6.9e-3,-strain_H(index_H),'r')
hold on
plot(Time_L(index_L)-tevent,strain_L(index_L),'g')
xlabel("time (s) since "+num2str(tevent))
ylabel('strain')
title('Advanced LIGO strain data near GW150914')
legend('H1 strain','L1 strain')
% fs = 4096;
% NFFT = fs;
% fmin = 10;
% fmax = 2000;
% [P_H,freqs_H] = pwelch(strain_H,fs);
% psd_H = interp1(freqs_H,P_H);
% loglog(freqs_H,sqrt(P_H))

%% 大作业
[time strain] = textread('observed-H.txt','%f %f');
[time_l,strain_l] = textread('observed-L.txt','%f %f','headerlines',1);
time = time-6.9e-3;
strain = -strain;
plot(time_l,strain_l,'LineWidth',1)
hold on
plot(time,strain,'LineWidth',1)
set(gca,'FontSize',12,'xTicklabel',{' ','0.3','0.35','0.4','0.45'})
set(gca,'yTicklabel',{' ','-1','-0.5','0','0.5','1','1.5'})
grid on
xlabel('Time (s)','FontSize',12)
ylabel('Strain(10^{-21})')
xlim([0.25 0.45])


%% 第一问
clc;
% clear;

f = @(theta,phi) abs(0.5*(1+cos(theta).^2).*cos(2*phi));
f_x = @(theta,phi) abs(cos(theta).*sin(2*phi));
theta = linspace(0,pi,200);
phi = linspace(0,2*pi,400);
[Theta,Phi] = meshgrid(theta,phi);
T = f(Theta,Phi);
[X,Y,Z] = ToRad(T,Theta,Phi);
mesh(X,Y,Z)
xlabel('x')
ylabel('y')
zlabel('z')
title('F_+(\theta,\phi)')
figure
T1 = f_x(Theta,Phi);
[X1,Y1,Z1] = ToRad(T1,Theta,Phi);
mesh(X1,Y1,Z1)
xlabel('x')
ylabel('y')
zlabel('z')
title('F_\times(\theta,\phi)')
T3 = sqrt(T1.^2+T.^2);
figure
[X3,Y3,Z3] = ToRad(T3,Theta,Phi);
mesh(X3,Y3,Z3)
xlabel('x')
ylabel('y')
zlabel('z')
title('$\sqrt{F^2_+(\theta,\phi)+F^2_\times(\theta,\phi)}$','interpreter','latex')
function [x,y,z]=ToRad(T,theta,phi)
x = T.*sin(theta).*cos(phi);
y = T.*sin(theta).*sin(phi);
z = T.*cos(theta);
end
