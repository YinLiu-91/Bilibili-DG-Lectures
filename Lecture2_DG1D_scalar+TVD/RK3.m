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
    
    % Stage I
    du = Lh(uh);
    uh1 = uh + dt*du;
    uh1 = TVD_Limiter(uh1);                     % TVD 每个euler步骤后都得加
    
    % Stage II
    du = Lh(uh1);
    uh2 = (3/4)*uh + (1/4)*uh1 + (1/4)*dt*du;
    uh2 = TVD_Limiter(uh2);                     % TVD 每个euler步骤后都得加
    
    % Stage III
    du = Lh(uh2);
    uh = (1/3)*uh + (2/3)*uh2 + (2/3)*dt*du;
    uh = TVD_Limiter(uh);                       % 经过TVD_Limiter才会称为真正的下一步的解的结果。TVD(Total Variation Diminishing，总变差递减)型限制器。
    
    fprintf('%d  %d\n',t,max(abs(uh(:,1))))
     
end