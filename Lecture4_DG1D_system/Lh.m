function du = Lh(uh)
global Nx NumGLP dimPk mm phiG phixG phiGR phiGL bcL bcR weight hx NumEq flux_type

% uhb是系数，不是真解
uhb = zeros(Nx + 2,dimPk,NumEq);
uhb(2:end - 1,:,:) = uh;

% uhG,uhR,uhL 都是真解，不是系数
uhG = zeros(Nx,NumGLP,NumEq);
uhR = zeros(Nx + 1,NumEq);      % 由于是方程组，所以多了一维NumEq
uhL = zeros(Nx + 1,NumEq);


fhat = zeros(Nx + 1,NumEq);
du = zeros(Nx,dimPk,NumEq);

% set_bc
if bcL == 1
    uhb(1,:,:) = uh(end,:,:);
elseif bcL == 2
    uhb(1,:,:) = uh(1,:,:);
end

if bcR == 1
    uhb(end,:,:) = uh(1,:,:);
elseif bcR == 2
    uhb(end,:,:) = uh(end,:,:);
end

% Step 1: calculate the Integral in cell
for i = 1:Nx
    for d = 1:dimPk
        for n = 1:NumEq       % 多了一个NumEq的循环
            for i1 = 1:NumGLP
                uhG(i,i1,n) = uhG(i,i1,n) + uh(i,d,n)*phiG(i1,d);
            end
        end
    end
end

for i = 1:Nx
    for d = 2:dimPk
        for i1 = 1:NumGLP % 由于向量方程，循环了三次
            du(i,d,1) = du(i,d,1) + 0.5*weight(i1)*f1(uhG(i,i1,1),uhG(i,i1,2),uhG(i,i1,3))*phixG(i1,d);
            du(i,d,2) = du(i,d,2) + 0.5*weight(i1)*f2(uhG(i,i1,1),uhG(i,i1,2),uhG(i,i1,3))*phixG(i1,d);
            du(i,d,3) = du(i,d,3) + 0.5*weight(i1)*f3(uhG(i,i1,1),uhG(i,i1,2),uhG(i,i1,3))*phixG(i1,d);
        end
    end
end

% Step 2: calculate the flux at edge
for i = 1:Nx + 1
    for d = 1:dimPk
        for n = 1:NumEq % 多了一个NumEq的循环
            uhR(i,n) = uhR(i,n) + uhb(i,d,n)*phiGR(d);
            uhL(i,n) = uhL(i,n) + uhb(i + 1,d,n)*phiGL(d);
        end
    end
end

for i = 1:Nx + 1
    uR = uhL(i,:);
    uL = uhR(i,:);
    % 下面都是根据通量的公式得来的，公式详见https://zhuanlan.zhihu.com/p/623709134
    fR(1) = f1(uR(1),uR(2),uR(3));
    fR(2) = f2(uR(1),uR(2),uR(3));
    fR(3) = f3(uR(1),uR(2),uR(3));
    fL(1) = f1(uL(1),uL(2),uL(3));
    fL(2) = f2(uL(1),uL(2),uL(3));
    fL(3) = f3(uL(1),uL(2),uL(3));
    [SLmax,SLmin] = wavespeed(uL);
    [SRmax,SRmin] = wavespeed(uR);
    SR = max(SLmax,SRmax);
    SL = min(SLmin,SRmin);
    switch flux_type
        case 1
            fhat(i,:) = LF_Flux(uR,uL,fR,fL,SR,SL);
        case 2
            fhat(i,:) = HLL_Flux(uR,uL,fR,fL,SR,SL);
        case 3
            fhat(i,:) = HLLC_Flux(uR,uL,fR,fL,SR,SL);
    end
end

% Step 3: calculate [f*phi] on each cell
for i = 1:Nx
    for d = 1:dimPk
        for n = 1:NumEq
            du(i,d,n) = du(i,d,n) - (1/hx)*(phiGR(d)*fhat(i + 1,n) - phiGL(d)*fhat(i,n));
        end
    end
end

for d = 1:dimPk
    du(:,d,:) = du(:,d,:)/mm(d);
end

end