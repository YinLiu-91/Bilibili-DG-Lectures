function du = Lh(uh)
global Nx NumGLP dimPk mm phiG phixG phiGR phiGL bcL bcR weight hx
uhb = [[0,0,0];uh;[0,0,0]];
uhG = zeros(Nx,NumGLP);
uhR = zeros(Nx + 1,1);
uhL = zeros(Nx + 1,1);
fhat = zeros(Nx + 1,1);
du = zeros(Nx,dimPk);

% set_bc 周期边界
if bcL == 1
    uhb(1,:) = uh(end,:);
end

if bcR == 1
    uhb(end,:) = uh(1,:);
end

% Step 1: calculate the Integral in cell
for i = 1:Nx              % 对所有对单元上的
    for d = 1:dimPk       % 所有型函数进行计算，并累加到uhG上。进行对是所有型函数对累加
        uhG(i,:) = uhG(i,:) + uh(i,d)*phiG(:,d)';
    end
end

for i = 1:Nx              % 对所有对单元上对
    for d = 2:dimPk       % 因为0阶为常数，导数为0 ，所以从2开始
        for i1 = 1:NumGLP % 所有积分点，计算积分结果
            du(i,d) = du(i,d) + 0.5*weight(i1)*f(uhG(i,i1))*phixG(i1,d);
        end
    end
end

% Step 2: calculate the flux at edge
for i = 1:Nx + 1          % 对所有对网格节点（不是有限元点，没有网格中间的那个点），注意这里的uhR，uhL是相对于单元来说的
    for d = 1:dimPk       % 计算所有的型函数
        uhR(i) = uhR(i) + uhb(i,d)*phiGR(d);          % 累加所有型函数的结果
        uhL(i) = uhL(i) + uhb(i + 1,d)*phiGL(d);      % 累加所有型函数的结果 而 flux 函数中的 u-，u+(程序中为 uL 和 uR )是相对于那个端点而言的左值和右值
    end
end

for i = 1:Nx + 1          % 循环所有节点
    uR = uhL(i);          % 计算uR，uR是相对于节点来说的，uhL(1)的起点是真实单元的左边值
    uL = uhR(i);          % uhR(1)是第一单元的上一个单元的右边值（即因周期边界虚拟的单元)，这里uhR比uhL往前一个记录位置感觉
    alpha = 1;            %
    fhat(i) = 0.5*(f(uR) + f(uL) - alpha*(uR - uL));  % 人工通量计算，从此看出，uR，uL是相对于节点来说的
end

% Step 3: calculate [f*phi] on each cell
for i = 1:Nx              % 对所有单元循环
    for d = 1:dimPk       % 累加所有基函数，下面看出，fhat(i+1)是单元的右端值，fhat(i)是单元的左端值
        du(i,d) = du(i,d) - (1/hx)*(phiGR(d)*fhat(i + 1) - phiGL(d)*fhat(i));
    end
end

for d = 1:dimPk         % 质量阵直接相除即可
    du(:,d) = du(:,d)/mm(d);
end

end
