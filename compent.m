function T = compent(C)
% 输入一组行向量，数据存在跳变，输出的向量是以第一个元素为基准补偿所有跳变的分量得到的最终值
l = length(C);
delta = 0;%补偿相位差
delta1 = 0; %补偿平衡量
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