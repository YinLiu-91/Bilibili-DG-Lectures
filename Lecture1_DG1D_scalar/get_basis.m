% get_basis.m
global phiG phixG mm phiGR phiGL
phiG = zeros(NumGLP,dimPk);         % 型函数,3个型函数为：1. \phi0=1;2. \phi1=x; 3. \phi2=x^2-1/3
phixG = zeros(NumGLP,dimPk);        % 型函数的导数
phiGR = zeros(1,dimPk);
phiGL = zeros(1,dimPk);
mm = zeros(1,dimPk);

for i = 1:NumGLP                    % 这里基函数在不同的积分点处不同的阶次基函数有不同的值(3*5=15)
    phiG(i,1) = 1;                  % 型函数点0次项值
    phiG(i,2) = lambda(i);          % 型函数点1次项值
    phiG(i,3) = lambda(i)^2 - 1/3;  % 型函数点2次项值

    phixG(i,1) = 0;                 % 型函数点0次项导数值
    phixG(i,2) = 1/hx1;             % 型函数点1次项导数值
    phixG(i,3) = 2*lambda(i)/hx1;   % 型函数点2次项导数值
end

% 下面计算了单元左右断点处的型函数的值，和系数乘后累加可用于计算左右断点处的解值
phiGR(1) = 1;                       % 型函数0在右端的值
phiGR(2) = 1;                       % 型函数1在右端的值
phiGR(3) = 2/3;                     % 型函数2在右端的值

phiGL(1) = 1;
phiGL(2) = -1;
phiGL(3) = 2/3;

mm(1) = 1;                          % 质量矩阵的(0,0)处值， 这里的值是对3个型函数在区间上的积分结果
mm(2) = 1/3;                        % 质量矩阵的(1,1)处值
mm(3) = 4/45;                       % 质量矩阵的(2,2)处值
