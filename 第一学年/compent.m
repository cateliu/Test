function T = compent(C)
% ����һ�������������ݴ������䣬������������Ե�һ��Ԫ��Ϊ��׼������������ķ����õ�������ֵ
l = length(C);
delta = 0;%������λ��
delta1 = 0; %����ƽ����
T(1) = C(1);
for j = 1:l-1
    chazhi = C(j+1)-C(j);
        if abs(chazhi)>=3
        if j == 1
            delta1 = C(4)-C(3);
%         elseif j==2
%             delta1 = C(4)-C(3);
        else
            delta1 = C(j)-C(j-1);
        end
        delta = delta+chazhi-delta1;
        end
    T(j+1) = C(j+1)-delta;
end
end