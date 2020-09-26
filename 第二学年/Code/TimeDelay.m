% % 探测相位为0-2pi
% % 时延最大为10ms
% a = 0:0.01:2*pi;
% x = 0:0.01:1;
% %计算第二帧的时延对相位影响

% [A,X] = meshgrid(a,x);
% C = sin(A+X/2).*cos(X/2)./cos(A);
% partial = (cos(A+X/2).*cos(X/2)-sin(X/2).*sin(A+X/2))./cos(A);
% y = 1/2*1./(1+C.^2).*partial;
% 
% 
% mesh(A,X,y)
% xlabel('A')
% ylabel('X')
% clear
% clc
% a = 0.5;
% syms a tau1 tau2 tau3 
% f(a, tau1, tau2, tau3 ) = atan(sin(a+(tau1+tau2)/2)*cos((tau3-tau1)/2)/(cos(a+tau2/2)*cos(tau2/2)));
% d1 = diff(f,tau1);
clear 
clc
C = @(a,tau1,tau2,tau3) (sin(a+(tau3+tau1)/2).*cos((tau3-tau1)/2))./(cos(a+tau2/2).*(cos((tau3-tau1)/2)));
y1 = @(a,tau1,tau2,tau3) 1./(1+C(a,tau1,tau2,tau3).^2);
partial1 = @(a,tau1,tau2,tau3) (cos(a+(tau3+tau1)/2).*cos((tau3-tau1)/2)-sin((tau1-tau3)/2).*sin(a+(tau3+tau1)/2))/(cos(a+tau2/2).*cos(tau2/2));
partial2 = @(a,tau1,tau2,tau3) (sin(a+(tau3+tau1)/2).*cos((tau3-tau1)/2).*(sin(a+tau2/2).*cos(tau/2)+sin(tau2/2).*cos(a+tau2/2)))./(cos(a+tau2/2).^2.*cos(tau2/2).^2);
P = @(a,tau1,tau2,tau3) 1/2*y1(a,tau1,tau2,tau3).*partial1(a,tau1,tau2,tau3);
P2 = @(a,tau1,tau2,tau3) 1/2*y1(a,tau1,tau2,tau3).*partial1(a,tau1,tau2,tau3);

tau = -1e-3:1e-4:1e-3;
a = 0:0.01:2*pi;

Max = P(a(1),tau(1),tau(1),tau(1));
Max2 = P2(a(1),tau(1),tau(1),tau(1));
%%
for j = 1:length(a)
    for tau1 = 1:length(tau)
        for tau2 = 1:length(tau)
            for tau3 = 1:length(tau)
                M(j,tau1,tau2,tau3) = P2(a(j),tau(tau1),tau(tau2),tau(tau3));
                if abs(M(j,tau1,tau2,tau3)) >= Max2
                    Max = M(j,tau1,tau2,tau3);
                    Xia = [j,tau1,tau2,tau3];
                end
            end
        end
    end
end
Xia
Max
[a(Xia(1)),tau(Xia(2)),tau(Xia(3)),tau(Xia(4))]
% 当延迟为-1s - 1s
% 相位 0 - 2pi
% Xia =315    21     1    21
% Max =0.5000
% ans =3.1400    0.0010   -0.0010    0.0010
%%
% m = max(M,[],'all');
% for j = 1:length(a)
%     for tau1 = 1:length(tau)
%         for tau2 = 1:length(tau)
%             for tau3 = 1:length(tau)
%                 if M(j,tau1,tau2,tau3) == m
%                     Xia = [ j,tau1,tau2,tau3];
%                     sprintf('%d %d %d %d', j,tau1,tau2,tau3);
%                 end
%             end
%         end
%     end
% end
