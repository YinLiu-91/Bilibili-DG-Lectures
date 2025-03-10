% get_basis.m
global phiG phiGLL phixG mm phiGR phiGL
phiG = zeros(NumGLP,dimPk);
phiGLL = zeros(NumGLP,dimPk);
phixG = zeros(NumGLP,dimPk);
phiGR = zeros(1,dimPk);
phiGL = zeros(1,dimPk);
mm = zeros(1,dimPk);

for i = 1:NumGLP
    phiG(i,1) = 1;                  % 这里基函数在不同的积分点处不同的阶次基函数有不同的值(3*5=15)
    phiG(i,2) = lambda(i);
    phiG(i,3) = lambda(i)^2 - 1/3;
    
    phiGLL(i,1) = 1;
    phiGLL(i,2) = lambdaL(i);
    phiGLL(i,3) = lambdaL(i)^2 - 1/3;
    
    phixG(i,1) = 0;
    phixG(i,2) = 1/hx1;
    phixG(i,3) = 2*lambda(i)/hx1;
end

phiGR(1) = 1;
phiGR(2) = 1;
phiGR(3) = 2/3;

phiGL(1) = 1;
phiGL(2) = -1;
phiGL(3) = 2/3;

mm(1) = 1;
mm(2) = 1/3;
mm(3) = 4/45;