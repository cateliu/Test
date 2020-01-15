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

plot(Time_H(index_H)-tevent,strain_H(index_H),'r')
hold on
plot(Time_L(index_L)-tevent,strain_L(index_L),'g')
xlabel("time (s) since "+num2str(tevent))
ylabel('strain')
title('Advanced LIGO strain data near GW150914')
legend('H1 strain','L1 strain')
%% 功率谱密度
figure
fs = 4096;
NFFT = 2048;
fmin = 10;
fmax = 2000;
[P_H,freqs_H] = pburg(strain_H,2047,0:2048,fs);
[P_L,freqs_L] = pburg(strain_L,2047,0:2048,fs);
% [P_L,freqs_L] = pwelch(strain_L,fs,NFFT,fs,fs);

psd_H =@(freqs) interp1(freqs_H,P_H,freqs);
psd_L =@(freqs) interp1(freqs_L,P_L,freqs);
loglog(freqs_H,sqrt(P_H),'r')
hold on
loglog(freqs_L,sqrt(P_L),'g')
ylabel('ASD(strain/rtHz)')
xlabel('Freq(Hz)')
legend('H1 strain','L1 strain')
title('Advanced LIGO strain data near GW150914')
axis([fmin fmax 1e-24 1e-19])
%%  白化 whitening
strain_H1_whiten = whiten(strain_H,psd_H,ts_H);
strain_L1_whiten = whiten(strain_L,psd_L,ts_L);
[bb,ab] = butter(4,[20*2./fs,300*2./fs],'bandpass');
strain_H1_whitenbp = filtfilt(bb,ab,strain_H1_whiten);

%% d带通滤波器
[bb,ab] = butter(4, [20./2./fs, 300./2./fs], 'bandpass');
y= bandpass(strain_H,[150,200], 4096);
plot(y)
axis([0 0.4*4096 -1.5e-21 1.5e-21])
strain_white = filtfilt(bb, ab, strain_H);

%% 大作业
clear
clc
close all
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

T_index = find(strain(2:end).*strain(1:end-1)<0);
T_index_t = T_index(find(time(T_index)>0.36&time(T_index)<0.426));
plot(time,strain,'.')
hold on
grid on
plot(time(T_index),strain(T_index),'--')
plot(time(T_index_t),strain(T_index_t),'+')
f = 1./2./(time(T_index_t(2:end))-time(T_index_t(1:end-1)));
f_1 = f.^(-8/3);
figure
plot(time(T_index_t(1:end-1)),f_1,'.')
ylabel('(Frequency/Hz)^{-8/3}')
xlabel('Time(s)')
title('Linear fit of f_{GW}^{-8/3}(t)')
%%
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
[E,Q] = meshgrid(log10(q),e);
contour(E,Q,R',0.5:0.25:2);
xlabel('Mass Ratio(q)')
ylabel('Eccentricity(e)')

