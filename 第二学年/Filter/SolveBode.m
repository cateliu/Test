subplot 211

semilogx(X__1(:,1), X__1(:,2))
xlim([1e-3,1e9])
title("Bodeͼ")
grid on
ylabel("����/dB")
subplot 212
semilogx(X__1(:,1), X__1(:,4))
xlim([1e-3,1e9])
ylim([-360 -90])
yticks([-360 -270 -180 -90])
ylabel("��λ/deg")
xlabel("Ƶ��/Hz")
grid on
