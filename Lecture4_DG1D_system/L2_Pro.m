% L2_Pro.m
uh = zeros(Nx,dimPk,NumEq);

for i = 1:Nx
    for d = 1:dimPk
        for n = 1:NumEq       % 这里相比之前只是多了一个对于变量的循环，其他没变
            for i1 = 1:NumGLP % https://zhuanlan.zhihu.com/p/605450187  
            uh(i,d,n) = uh(i,d,n) + 0.5*weight(i1)*ureal(i,i1,n)*phiG(i1,d); % 这里累加的是积分结果，对于基函数来说，在不同积分点处不同的阶数的基函数有不同的值
            end
        end
    end
end

for d = 1:dimPk
    uh(:,d,:) = uh(:,d,:)/mm(d);
end