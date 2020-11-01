%phi应该是一维数组，judge是判断相邻像元的是否跳变的临界值，一般设为3
%jump为相位跳变的值，通常为pi
function Phi = PhaseUnwrapping(phi,judge,jump)
%     x = 1:0.01:10;
%     phi = mod(x,jump);
    [m,n] = size(phi);
    if isvector(phi)
        plot(phi);
        v = vertical(phi, judge);
        c0 = 0;
        c = count(v,c0);
        figure
        % 判断行向量还是列向量
        if isrow(phi)
            Phi = phi+jump*c;
        else
            Phi = phi+jump*c';
        end
        plot(Phi)
        
    elseif m>1 && n>1
        Phi = TPhaseUnwrapping(phi,judge, jump);
        mesh(Phi)
    end
    
end
function Phi = TPhaseUnwrapping(phi,judge, jump)
    [m,n] = size(phi);

    v = vertical(phi(:,1), judge);

    for i = 1:m
        z(i,:) = vertical(phi(i,:), judge);
    end
%     c = zeros(m,n);
    c(:,1) = count(v,0);
    for j = 2:n
        c(:,j) = z(:,j-1) + c(:,j-1);
    end
    
    Phi = phi + jump*c;
end

function v = vertical(phi, judge)
    for j = 1:length(phi)-1
        delta = phi(j+1)-phi(j);
        if delta >0 && floor(abs(delta)/judge)
            v(j) = -1;
        elseif delta <0 && floor(abs(delta)/judge)
            v(j) = 1;
        else
            v(j) = 0;
        end
    end
end

function c = count(v,c0)
    for j = 1:length(v)+1
        if j == 1
            c(j) = c0;
        else
            c(j) = c(j-1) + v(j-1);
        end
    end
end
