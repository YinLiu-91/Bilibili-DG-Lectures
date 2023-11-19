% RK3.m

dt = CFL*hx;
t = 0;
% add for plot
nsaveT = 2
saveT = dt*nsaveT
saved=zeros(floor(tend/saveT)+1,Nx,dimPk);
i=0;

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
    
    %fprintf('%d  %d\n',t,max(abs(uh(:,1))))
    
    % add for plot
    if(mod(floor(t/dt),nsaveT)==0)
      i=i+1;
      % 在单元[-1,0,1]的三个位置处计算最终的求解变量结果
      saved(i,:,:)=[uh(:,1)-uh(:,2)+2/3*uh(:,3) uh(:,1)-1/3*uh(:,3) uh(:,1)+uh(:,2)+2/3*uh(:,3)];
    end
     
end

% plot 
% get coordinate
% 
end_coor=xb;
xNode=0:end_coor/Nx:end_coor;       % 生成单元边界点
xNode1=xNode+end_coor/Nx/2;  % 生成单元中心点
%xNode1=1/Nx:1/Nx;
xh = reshape(sort([xNode,xNode(2:end-1),xNode1(1:end-1)]),dimPk,Nx); % 这里是为了产生单元的左中右点坐标位置，因为单元左右端点在不同单元有不同值，所以这样设置坐标
h=figure;set(gcf, 'Color','white')
for i=1:length(saved)
    y=reshape(saved(i,:,:),3,100);
    plot(xh,y)
    axis([0 end_coor -1.5 1.5])
    pause(.01)
end