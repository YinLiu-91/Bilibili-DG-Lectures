% RK3.m

dt = CFL*hx;
t = 0;

while t < tend

    if t + dt >= tend
        dt = tend - t;
        t = tend;
    else
        t = t + dt;
    end

    % Stage I  前向欧拉；uh刚开始是我们给对初始值
    du = Lh(uh);
    uh1 = uh + dt*du;

    % Stage II 对中间值再使用前向欧拉
    du = Lh(uh1);
    uh2 = (3/4)*uh + (1/4)*uh1 + (1/4)*dt*du;

    % Stage III 对中中间值再再使用前向欧拉
    du = Lh(uh2);
    uh = (1/3)*uh + (2/3)*uh2 + (2/3)*dt*du;

    fprintf('%d  %d\n',t,max(abs(uh(:,1))))

end
