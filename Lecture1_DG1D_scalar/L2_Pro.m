% L2_Pro.m
% 为什么要用L2投影获取初值，这里回答了：https://zhuanlan.zhihu.com/p/604524727
% 因为要得到0时刻此时的u在每个单元的基函数的系数c0，然后再对c0进行时间推进。
% 本质上是通过某一时刻u在每个单元的基函数的系数，刻画这个时刻的计算域内u的的分布。

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
