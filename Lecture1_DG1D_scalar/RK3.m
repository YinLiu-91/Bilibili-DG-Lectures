% RK3.m

dt = CFL*hx;
t = 0;
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

    % Stage I  前向欧拉；uh刚开始是我们给对初始值
    du = Lh(uh);
    uh1 = uh + dt*du;

    % Stage II 对中间值再使用前向欧拉
    du = Lh(uh1);
    uh2 = (3/4)*uh + (1/4)*uh1 + (1/4)*dt*du;

    % Stage III 对中中间值再再使用前向欧拉
    du = Lh(uh2);
    uh = (1/3)*uh + (2/3)*uh2 + (2/3)*dt*du;

    %fprintf('%d  %d\n',t,max(abs(uh(:,1))));
    if(mod(floor(t/dt),nsaveT)==0)
      i=i+1;
      % 在单元[-1,0,1]的三个位置处计算最终的求解变量结果
      %saved(i,:,:)=[uh(:,1)-uh(:,2)+2/3*uh(:,3) uh(:,1)-1/3*uh(:,3) uh(:,1)+uh(:,2)+2/3*uh(:,3)];
      saved(i,:,:)=uh;
    end

end

% plot 
% get coordinate
% 
end_coor=xb;
xNode=0:end_coor/Nx:end_coor;       % 生成单元边界点
xNode1=xNode+end_coor/Nx/2;  % 生成单元中心点
%xNode1=1/Nx:1/Nx;
xh = reshape(sort([xNode,xNode(2:end-1),xNode1(1:end-1)]),dimPk,Nx);
h=figure;set(gcf, 'Color','white')
for i=1:length(saved)
    y=reshape(saved(i,:,:),3,100);
    plot(xh,y)
    axis([0 end_coor -1.5 1.5])
%    if (i-1)/NormFreq == floor((i-1)/NormFreq)
%        j=j+1;
%    end
%    text(1.02,0.1,'L2-Norm')
%    text(1.02,0,num2str(norm2(j)));
%    text(1.02,-0.2,'RMS')
%    #text(1.02,-.3,num2str(rms(reshape(saved(2,:,i),K,1))));
%    text(1.02,1.1,'Time')
%    text(1.02,1,num2str((i-1)*saveT));
    %writeVideo(vidObj, getframe(h));
    pause(.01)
end
