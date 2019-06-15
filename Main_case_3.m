if 0
    getequations_model_1
    model_1
end
global zlaminitial;
global g gz gzd gu gud taoz taouz
global f fx fxd fu fud fr tao taou theta;
global impfun impfz impfzd impfu impfud taoimp taouimp;
global dimtheta dimtao dimx dimu dimrd;
global measure_system;
global xhistory
global uhistory
model_id=2;
switch model_id;
    case 1
        load('model_1.mat');
        xhistory=[2,0]';
        uhistory=[1;1];
        theta=[1,1,1,1];tao=1;taou=0;parameter=[theta,tao,taou];
        hx=[1 0;0 1];
        invQ=diag((1/0.3./[ 6.542740624338593 14.111654509129362]).^2);%检测噪声
        uminvalue=[0.05;5];
        umaxvalue=[0.20;35];xmin=[1,0];xmax=[10,0.0000001];
        rdvalue=zeros(dimrd,1);
    case 2
        load('model_2.mat');
        xhistory=[3.3e-4,4e-3]';
        uhistory=[1e+5;4.4e+5];
        theta=[1,1];tao=1;taou=0;parameter=[theta,tao,taou];
        hx=[1 0;0 1];
        invQ=diag((1/0.3./[ 3.3e-4 4e-3]).^2);%检测噪声s
        uminvalue=[0.6e+5;4.1e+5];
        umaxvalue=[1.5e+5;4.7e+5];
        rdvalue=[6e-4 9e-3]';
end
%%       
T=8;
dimy=size(hx,1);
Nsample=32;Nu=16;
tsample=linspace(0,T,Nsample+1);tsample=tsample(2:end);%采样点
measure_system=struct('tsample',tsample,'ysample',zeros(dimy,numel(tsample)),'hx',hx,'invQ',invQ);
tseq=linspace(0,T,Nu+1); 
% tt=linspace(0,T,200);dt=4*T/numel(tt);tt1=[];for i=1:numel(tseq);tt1=[tt1,linspace(max(0,tseq(i)-dt),min(T,tseq(i)+dt),10)];end;tt=unique(union(tt,tt1(tt1>=0&tt1<=T)));
tseq_rd=[0,tsample];rdseq_ref=rdvalue*ones(1,numel(tseq_rd)-1);rdfun_ref=@(t)interp_piecewise(tseq_rd,rdseq_ref,t);%参考干扰信号
rdseq=rdseq_ref;rdfun=@(t)interp_piecewise(tseq_rd,rdseq,t);
pestimate=parameter*1.0;pestimate(dimtheta+dimtao+1:end)=parameter(dimtheta+dimtao+1:end);%这个是已知的最好辨识结果，用于过程优化和信息优化
umin=uminvalue*ones(1,Nu);%umin=umin(:);
umax=umaxvalue*ones(1,Nu);%umax=umax(:);
%%
global len_p;%定义辨识实验设计中考虑的待辨识参数的长度
global criterion;
criterion='D';mode=3;%有路径约束优化模式--对于模型4，A判据很不凸，而D判据可以
len_p=dimtheta+dimtao;
%仿真
ubinary=umin;ubinary(1,2:2:end)=umax(1,2:2:end);ubinary(2,1:2:end)=umax(2,1:2:end);
ubk=umin+(umax-umin).*(randi(2,2,Nu)-1);
tic,[~,~,~,li,H]=ObjForInfoOptim(ubk,tseq,rdfun,parameter,xhistory,uhistory);toc
for i=1:dimx;
    figure;plot(li.t(li.t>=0),li.x(i,li.t>=0));
    prctile(li.x(i,:),80)
end
%优化1
len_p=dimtheta+dimtao;s=len_p;MaxIter=500;
tic,[uoptinfo1,fopt,u_kexing,xopt_b]=opt_Information_con(pestimate,rdfun,tseq,umin(:),umax(:),xhistory,uhistory,mode,MaxIter); toc
uopt_b=reshape(uoptinfo,2,Nu);
tic,[~,~,~,li_b,H]=ObjForInfoOptim(uopt_b,tseq,rdfun,parameter,xhistory,uhistory);toc
for i=1:dimx;
    figure;plot(li_b.t(li_b.t>=0),li_b.x(i,li_b.t>=0));
    prctile(li_b.x(i,:),80)
end

%优化2
len_p=dimtheta;s=len_p;
tic,[uoptinfo,fopt,u_kexing,xopt_a]=opt_Information_con(pestimate,rdfun,tseq,umin(:),umax(:),xhistory,uhistory,mode,MaxIter); toc
uopt_a=reshape(uoptinfo,2,Nu);
save('Main_case_3.mat');
%%
%检验辨识精度
ubinary=umin;ubinary(1,2:2:end)=umax(1,2:2:end);ubinary(2,1:2:end)=umax(2,1:2:end);
udict={ubinary,uopt_a,uopt_b};
xdict={xhistory,xhistory,xhistory};
cpdict={};
for kss=1:3;
    x=xdict{kss};
    u=udict{kss};
    [~,yout]=simulation(u,tseq,rdseq,tseq_rd,parameter,x,uhistory);
    cpopt=[];
    maxk=2;
    for k=1:maxk;
        noiseVar=sqrt(diag(inv(invQ)));%只考虑不相关的信号
        dimy=size(hx,1);ysample=yout;
        for i=1:dimy;ysample(i,:)=yout(i,:)+max(min(3*noiseVar(i),noiseVar(i)*randn(1,size(yout,2))),-3*noiseVar(i));end;
        popt=identify(ysample,u,rdseq_ref,s,pmin,pmax,pestimate,tseq,tseq_rd,xhistory,uhistory);
        cpopt=[cpopt;popt(:)'];
    end;
    cpdict{kss}=cpopt;
end

%%
%%
%统计
for i=1:3;
    cpopt=cpdict{i};
    stdi=std(cpopt-ones(maxk,1)*ptrue(:)',1)%统计误差数据
end
[~,~,~,li,H]=ObjForInfoOptim(ubinary,tseq,rdfun,parameter,xhistory,uhistory);
[~,~,~,liopt_a,Hopt_a]=ObjForInfoOptim(uopt_a,tseq,rdfun,parameter,xopt_a,uhistory);
[~,~,~,liopt_b,Hopt_b]=ObjForInfoOptim(uopt_b,tseq,rdfun,parameter,xopt_b,uhistory);
fa=@(H)trace(inv(H(1:dimtheta+dimtao,1:dimtheta+dimtao)));
fe=@(H)sqrt(eig(inv(H(1:dimtheta+dimtao,1:dimtheta+dimtao))));
[fa(H),fa(Hopt_a),fa(Hopt_b)]
[fe(H),fe(Hopt_a),fe(Hopt_b)]
%作图
figure;
for i=1:2
    subplot(2,1,i);hold on;
    plot(li.t(li.t>=0),li.x(i,li.t>=0));
    plot(liopt_a.t(liopt_a.t>=0),liopt_a.x(i,liopt_a.t>=0),'r-')
    plot(liopt_b.t(liopt_b.t>=0),liopt_b.x(i,liopt_b.t>=0),'k--')
    if i==1;ylabel('x_{1}(m^2)');
    else;ylabel('x_{2}(m^2)');
    end
    xlabel('Time(h)');legend('BS','AD','AN'); 
end
endfigure;
subplot(1,2,1);plot_piecewise(tseq,ubinary(1,:)','r.');ylabel('u[m2]');
axis([0,T,umin(1),umax(1)*1.1])
subplot(1,2,2);plot_piecewise(tseq,ubinary(2,:)','r.');ylabel('v[m2]');
axis([0,T,umin(2),umax(2)*1.1])
title('Binary Input');
figure;
subplot(1,2,1);plot_piecewise(tseq,uopt_a(1,:)','r.');ylabel('u[m2]');
axis([0,T,umin(1),umax(1)*1.1])
subplot(1,2,2);plot_piecewise(tseq,uopt_a(2,:)','r.');ylabel('v[m2]');
axis([0,T,umin(2),umax(2)*1.1])
title('A-optimal Input (without delay)')
figure;
subplot(1,2,1);plot_piecewise(tseq,uopt_b(1,:)','r.');ylabel('u[m2]');
axis([0,T,umin(1),umax(1)*1.1])
subplot(1,2,2);plot_piecewise(tseq,uopt_b(2,:)','r.');ylabel('v[m2]');
axis([0,T,umin(2),umax(2)*1.1])
title('A-optimal Input (with delay)')

udict={ubinary,uopt_a,uopt_b};

for kss=1:3;
    u=udict{kss};
    [~,yout]=simulation(u,tseq,rdseq,tseq_rd,parameter,xhistory,uhistory);
    cpopt=[];
    for k=1:20;
        noiseVar=sqrt(diag(inv(invQ)));%只考虑不相关的信号
        dimy=size(hx,1);ysample=yout;
        for i=1:dimy;ysample(i,:)=yout(i,:)+max(min(3*noiseVar(i),noiseVar(i)*randn(1,size(yout,2))),-3*noiseVar(i));end;
        popt=identify(ysample,u,rdseq_ref,s,pmin,pmax,pestimate,tseq,tseq_rd,xhistory,uhistory);
        cpopt=[cpopt;popt(:)'];
    end;
    cpdict{kss}=cpopt;
end
% cpopt_A=cpopt;
% save('cpopt_A.mat','cpopt_A')
