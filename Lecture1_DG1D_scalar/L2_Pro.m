% L2_Pro.m
uh = zeros(Nx,dimPk);         % [100x3]大小，因为是3阶次，所以单元中点也有一个系数

for i = 1:Nx
    for d = 1:dimPk
        for i1 = 1:NumGLP     % 对所有对G-L积分点进行累加，\int {u*\phi}
            uh(i,d) = uh(i,d) + 0.5*weight(i1)*ureal(i,i1)*phiG(i1,d);  % phiG存储了3个型函数在积分点处点值
        end
    end
end
% plot(uh(:,1))
for d = 1:dimPk
    uh(:,d) = uh(:,d)/mm(d);  % mm在get_basis.m中定义, 这里看出，uh只是多项式的系数，并不是解
end
