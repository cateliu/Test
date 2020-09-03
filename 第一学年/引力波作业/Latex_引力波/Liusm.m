clear
clc
close all 
[time strain] = textread('observed-H.txt','%f %f');
[time_l,strain_l] = textread('observed-L.txt','%f %f','headerlines',1);
subplot 221
plot(time_l,strain_l,'LineWidth',1)
hold on
plot(time,strain,'LineWidth',1)%图1
set(gca,'FontSize',12,'xTicklabel',{' ','0.3','0.35','0.4','0.45'})
set(gca,'yTicklabel',{' ','-1','-0.5','0','0.5','1','1.5'})
grid on
xlabel('Time (s)','FontSize',12)
ylabel('Strain(10^{-21})')
title('未修正的引力波波形')
xlim([0.25 0.45])

time = time-6.9e-3;
strain = -strain;
subplot 222
plot(time_l,strain_l,'LineWidth',1)
hold on
plot(time,strain,'LineWidth',1)%图1
set(gca,'FontSize',12,'xTicklabel',{' ','0.3','0.35','0.4','0.45'})
set(gca,'yTicklabel',{' ','-1','-0.5','0','0.5','1','1.5'})
grid on
xlabel('Time (s)','FontSize',12)
ylabel('Strain(10^{-21})')
title('修正的引力波波形')
xlim([0.25 0.45])

T_index = find(strain(2:end).*strain(1:end-1)<0);
T_index_t = T_index(find(time(T_index)>0.36&time(T_index)<0.426));
f = 1./2./(time(T_index_t(2:end))-time(T_index_t(1:end-1)));
f_1 = f.^(-8/3);
subplot 223
plot(time(T_index_t(1:end-1)),f_1,'.')%图3
hold on
a = polyfit(time(T_index_t(1:end-1)),f_1,1);
plot(time(T_index_t(1:end-1)),polyval(a,time(T_index_t(1:end-1))))
ylim([0 4e-5])
ylabel('(Frequency/Hz)^{-8/3}')
xlabel('Time(s)')
title('Linear fit of f_{GW}^{-8/3}(t)')
%% 图5
e = 0:0.001:0.8;
q = 1:0.1:100;
[E,Q] = meshgrid(q,e);
l = (1-e.^2).^(-7/2).*(1+73/24*e.^2+37/96*e.^4);
M0 = 2e30;
c = 3e8;
G = 6.67e-11;
fmax = 150;
M = l.^(-3/5)*M0;
q1 = q';
R = (1-e).*l.^(2/5).*2.66.*q1.^(2/5)./(1+q1).^(4/5);
subplot 224
[E,Q] = meshgrid(log10(q),e);
contour(E,Q,R',0.5:0.25:2);
xlabel('Mass Ratio(q)')
ylabel('Eccentricity(e)')
