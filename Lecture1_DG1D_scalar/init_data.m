% init_data.m
global bcL bcR hx hx1

%----------------
xa = 0;             %  求解域起点
xb = 2*pi;          %  求解域终点
u0 = @(x) sin(x);
bcL = 1;
bcR = 1;
tend = 2*pi;
%----------------

hx = (xb - xa)/Nx;          % 单元的长度
hx1 = 0.5*hx;               % 单元长度的一半

Xc = xa + hx1:hx:xb - hx1;  % 单元的中点坐标

ureal = zeros(Nx,NumGLP);   % 最重要的是这个 ureal，它是用来储存初值(也是真解)在每个积分点上的值的。对每个固定的 i
                            % ureal(i,:)就是真解在第 i 个单元上的所有积分点上的值。
for i = 1:Nx
    for j = 1:NumGLP
        x=Xc(i) + hx1*lambda(j)
        ureal(i,j) = u0(Xc(i) + hx1*lambda(j)); % lambda 见get_GLP.m 是位置坐标
    end
end
