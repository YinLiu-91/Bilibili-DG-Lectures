% RK3.m

alpha = 1;
for i = 1:Nx
    [alpha1,~] = wavespeed(uh(i,1,:));
    if alpha1 > alpha
        alpha = alpha1;
    end
end

dt = CFL*hx/alpha;
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
    if Limit_type == 1
        uh1 = TVD_Limiter_P1(uh1);
    elseif Limit_type == 2
        uh1 = TVD_Limiter_P2(uh1);
    end
    %uh1 = pp_Limiter(uh1);
    
    % Stage II
    du = Lh(uh1);
    uh2 = (3/4)*uh + (1/4)*uh1 + (1/4)*dt*du;
    if Limit_type == 1
        uh2 = TVD_Limiter_P1(uh2);
    elseif Limit_type == 2
        uh2 = TVD_Limiter_P2(uh2);
    end
    %uh2 = pp_Limiter(uh2);
    
    % Stage III
    du = Lh(uh2);
    uh = (1/3)*uh + (2/3)*uh2 + (2/3)*dt*du;
    if Limit_type == 1
        uh = TVD_Limiter_P1(uh);
    elseif Limit_type == 2
        uh = TVD_Limiter_P2(uh);
    end
    %uh = pp_Limiter(uh);
    
    fprintf('%d  %d  %d\n',t,min(uh(:,1,1)),max(abs(uh(:,1,1))))
     
end