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
[P_H,freqs_H] = pwelch(strain_H,fs,NFFT,fs,fs);
[P_L,freqs_L] = pwelch(strain_L,fs,NFFT,fs,fs);
freqs_H1 = min(freqs_H):(max(freqs_H)-min(freqs_H))/length(min(freqs_H))/2:max(freqs_H);
freqs_L1 = min(freqs_L):(max(freqs_L)-min(freqs_L))/length(min(freqs_L))/2:max(freqs_L);

psd_H = interp1(freqs_H,P_H,freqs_H1);
psd_L = interp1(freqs_L,P_L,freqs_L1);
loglog(freqs_H,sqrt(P_H),'r')
hold on
loglog(freqs_L,sqrt(P_L),'g')
ylabel('ASD(strain/rtHz)')
xlabel('Freq(Hz)')
legend('H1 strain','L1 strain')
title('Advanced LIGO strain data near GW150914')
axis([fmin fmax 1e-24 1e-19])
%%  白化 whitening
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
figure
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
%%


figure
fs = 4096;
window = hamming(fs/16);
% [S F T P] = ;
strain_w = bandpass(strain,[35,350], 4096);
[S F T P] = spectrogram(strain_w,4096/16,4096/16-1,fs/16,fs,'yaxis');
mesh(T,F,P)
axis([0,0.4 0 350])
function strain_whiten=whiten(strain, interp_psd, dt)
Nt = length(strain);


end
