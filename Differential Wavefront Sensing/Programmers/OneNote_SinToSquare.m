clc
clear
close
% 画出SED上的电压信号强度曲线-正弦曲线
t = 0:0.01:30;
Sin = 3 + sin(t);
plot(t,Sin,'k', 'LineWidth',2);  %k表示黑色

% 画出经过比较器(Comparator)后形成的方波
Square = Sin;
Square(find(Sin>=3)) = 1;
Square(find(Sin<3)) =-1;
hold on
plot(t, Square, 'k', 'LineWidth',2)

% 方波信号经过FPGA后会产生触发信号,脉冲宽度为0.5
% 二分频
Dif = [0,diff(Sin)];
y1 = Sin;
y1(find(Dif >= 0 & Sin >= 3)) = -2;
y1(find(Dif >= 0 & Sin < 3)) = -4;
y1(find(Dif < 0 & Sin >= 3)) = -4;
y1(find(Dif < 0 & Sin < 3)) = -2;
hold on
plot(t, y1, 'k','LineWidth',2)

% 过零检测
temp1 = (y1(1:end-1)+3).*(y1(2:end)+3);
pulse_t = find(temp1 < 0);

% 脉冲触发信号
pulse = 0*t;
for i = 1:length(t)
    if ismember(i, pulse_t) 
        for j = 0:49
           pulse(i+j) = 1; 
        end
    end
end
hold on
plot(t, pulse(1:length(t))-6, 'k','LineWidth',2)

% 画图，对准采样点
y2 = -6:0.1:4;
[X,Y] = meshgrid(y2, t(pulse_t));
hold on
plot(Y',X','k--')
% 虚线
hold on
Virt_line1 = 0*t + 3;
Virt_line2 = 0*t + 0;
Virt_line3 = 0*t - 3;
plot(t, Virt_line1, 'k--');
plot(t, Virt_line2, 'k--');
plot(t, Virt_line3, 'k--');

ylim([-7,5]);
xticks([]);
yticks([]);
