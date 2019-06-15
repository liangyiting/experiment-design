function [uopt,fopt,u_kexing,xopt]=opt_Information(parameter,rdfun,tseq,umin,umax,xhistory,uhistory,mode,MaxIter)
xopt=xhistory;
options=optimset('display','iter','algorithm','interior-point','GradObj','on','gradconstr','on','TolFun',1e-14,'TolCon',1e-14,'TolX',1e-6);
% options.PlotFcns={@optimplotfunccount,@optimplotfval,...
% @optimplotfirstorderopt,@optimplotx};%
options.PlotFcns={@optimplotfval,@optimplotfirstorderopt,@optimplotfunccount,@optimplotconstrviolation};%@optimplotconstrviolation
try options.MaxIter=MaxIter;
catch options.MaxIter=500;
end
if mode==3 %有路径约束优化
    FUN=@(u)infoOptim(u,tseq,rdfun,parameter,xhistory,uhistory);
    tt1=linspace(tseq(1),tseq(end),floor(160*tseq(end)));%采用足够细的网格可以提高计算约束条件精度
    cfun=@(u)cfun_(u,tseq,tt1,rdfun,parameter,xhistory,uhistory);
    u0=umin+rand(size(umax)).*(umax-umin);if 0;[a,b]=check_diff(FUN,u0,1,1e-6);end
    if 0;f=@(u)newfun(u,tseq,tt1,rdfun,parameter,xhistory,uhistory);[a,b]=check_dconu(f,u0,1e-6);end%检测约束对u的梯度
    u_kexing=[];
    [uopt,fopt]=fmincon(FUN,umax,[],[],[],[],umin,umax,cfun,options);
elseif mode==4 %无路径约束，但是可对初始值进行优化
    FUN=@(u)infoOptim_x0(u,tseq,rdfun,parameter,uhistory);
    u0=umin+rand(size(umax)).*(umax-umin);
    if 0;[gdiff,df]=check_diff(FUN,u0,1,ep);end
    if 0; [gdiff,df]=check_dJx0(reshape(u0,2,numel(u0)/2),tseq,rdfun,parameter);end%检测目标值对初始状态的梯度是否准确
    ux0=[xhistory(:);u0(:)]';
    %xmin=0*xhistory;xmax=10*xhistory;
    xmin=[1,0];xmax=[10,0.0000001];%模型1的初始点范围
    uxmin=[xmin(:);umin(:)]';uxmax=[xmax(:);umax(:)]';
    [uxopt,fopt]=fmincon(FUN,ux0,[],[],[],[],uxmin,uxmax,[],options);
    u_kexing=[];
    global dimx;
    xopt=uxopt(1:dimx);
    uopt=uxopt(dimx+1:end);
elseif 1 %无路径约束，仅仅优化控制序列
    FUN=@(u)infoOptim(u,tseq,rdfun,parameter,xhistory,uhistory);
    tt1=linspace(tseq(1),tseq(end),floor(160*tseq(end)));%采用足够细的网格可以提高计算约束条件精度
    u0=umin+rand(size(umax)).*(umax-umin);if 0;[a,b]=check_diff(FUN,u0,1,1e-6);end
    u_kexing=[];
    [uopt,fopt]=fmincon(FUN,u0,[],[],[],[],umin,umax,[],options);
end
function [J,dJ]=kexingfun(u)
J=0;dJ=0*u;
function [gdiff,df]=check_dconu(fun,useq,ep)
%检测约束对u的梯度是否准确
[f0,df]=fun(useq);dimu=size(useq,1);gdiff=0*df;
try ep;catch ep=1e-6;end
for i=1:dimu%numel(useq);
    ui=useq;ui(i)=useq(i)+ep;gdiff(i,:)=(fun(ui(:))'-f0')/ep;
end;
format long
try;
    [gdiff(1:i,1),df(1:i,1)]
end;
function [gdiff,df]=check_dJx0(useq,tseq,rdfun,parameter)
global xhistory;
fun=@(xhistory)ObjForInfoOptim(useq,tseq,rdfun,parameter,xhistory);
x=xhistory;x=x.*rand(size(x));
[f0,~,df]=fun(x);
gdiff=0*df;
ep=1e-6;
for i=1:numel(x);
    xi=x;xi(i)=x(i)+ep;gdiff(i)=(fun(xi)'-f0')/ep;
end;
format long
[gdiff,df]
function [c,gc]=newfun(u,tseq,tt,rdfun,parameter,xhistory,uhistory);
[c,ceq,gc,gceq]=cfun_(u,tseq,tt,rdfun,parameter,xhistory,uhistory);
function [c,ceq,gc,gceq]=cfun_(u,tseq,tt,rdfun,parameter,xhistory,uhistory)
global dimu;
T=tseq(end);useq=reshape(u,dimu,numel(u)/dimu);
ceq=[];gceq=[];
c=[];gc=[];
% epsilon1=0.00001;epsilon2=0.00003;
if 0
    x1max=0.00090;epsilon1=0.000001;
    x2max=0.0070;epsilon2=0.00001;
else
    global x1max x2max epsilon1 epsilon2;
%     x1max=30;epsilon1=0.1;
%     x2max=10;epsilon2=0.1;
end
consystem1=struct('L',@(t,x,u)phi(x(1)-x1max,epsilon1),'Lx',@(t,x,u)dphi(x(1)-x1max,epsilon1)*[1;0],'Lu',@(t,x,u)[0,0]');
if nargout>3;
    [c1,gc1,tt,xx]=ObjForProcessOptim(useq,tseq,tt,rdfun,parameter,consystem1,xhistory,uhistory);
else
    c1=ObjForProcessOptim(useq,tseq,tt,rdfun,parameter,consystem1,xhistory,uhistory);
end
consystem2=struct('L',@(t,x,u)phi(x(2)-x2max,epsilon2),'Lx',@(t,x,u)dphi(x(2)-x2max,epsilon2)*[0;1],'Lu',@(t,x,u)[0,0]');
if nargout>3;
    [c2,gc2,tt,xx]=ObjForProcessOptim(useq,tseq,tt,rdfun,parameter,consystem2,xhistory,uhistory);
else
    c2=ObjForProcessOptim(useq,tseq,tt,rdfun,parameter,consystem2,xhistory,uhistory);
end
% c=[c;c1-epsilon1*T-0.002;c2-epsilon2*T-0.019]';%路径约束
c=[c;c1-epsilon1*T*0-0.000;c2-epsilon2*T*0-0.00]';%路径约束
if nargout>3;gc=[gc(:)';gc1(:)';gc2(:)']';end
% if 0; li=[1,0.5];for i=1:dimu;[ci,gci]=smoothFun(useq,i);c=[c,ci/1e10-li(i)];gc=[gc,gci(:)/1e10];endend
K=1e0;
c=c(:)*K;
if nargout>2; gc=gc*K;end
function cnew=phi(c,epsilon);
if c>epsilon;cnew=c;
elseif c<-epsilon;cnew=0;
else cnew=(c+epsilon).^2/4/epsilon;
end;
function dc=dphi(c,epsilon);
if c>epsilon; dc=1;
elseif c<-epsilon;dc=0;
else dc=(c+epsilon)/2/epsilon;
end;
function [c3,gc3]=smoothFun(useq,i)
dimu=size(useq,1);
li=diff(useq(i,:));
c3=1/2*norm(li)^2;
li1=[0,li];
li2=[li,0];li3=li1-li2;
gc3=zeros(dimu,size(useq,2));gc3(i,:)=li3;
function [J,dJ,tt,xx,ttlam,clam,grad]=ObjForProcessOptim(useq,tseq,tt,rdfun,parameter,pro_system,xhistory,uhistory)
global g gx gxd gu gud taoz taouz
global f fx fxd fu fud fr;
global impfun impfz impfzd impfu impfud taoimp taouimp;
global dimtheta dimtao dimx dimu;
tao=parameter(dimtheta+1:dimtheta+dimtao);taou=parameter(dimtheta+dimtao+1:end);
history.x=[xhistory(:),xhistory(:)];
history.t=[-20,0];
history.u=uhistory(:)*ones(numel(history.t)-1);%初始的控制应该是已知的
state_system_u=struct('f',@(t,x,xd,u,ud)f(t,x,xd,u,ud,rdfun(t),parameter),...
    'fx',@(t,x,xd,u,ud)fx(t,x,xd,u,ud,rdfun(t),parameter),...
    'fxd',@(t,x,xd,u,ud)fxd(t,x,xd,u,ud,rdfun(t),parameter),...
    'fu',@(t,x,xd,u,ud)fu(t,x,xd,u,ud,rdfun(t),parameter),...
    'fud',@(t,x,xd,u,ud)fud(t,x,xd,u,ud,rdfun(t),parameter),...
    'tao',tao,'taou',taou,'dimx',dimx);
useq1=[history.u,useq];
tseq1=unique([history.t,tseq]);

if nargout>1;
    if isstruct(pro_system)
        [J,dJ,tt,xx,ttlam,clam,grad]=idenObj(useq1,tseq1,tt,history,state_system_u,pro_system,[]);
        dJ=dJ(:,2:end);
    else
        [~,~,tt,xx]=idenObj(useq1,tseq1,tt,history,state_system_u,[],[]);
        J=[];dJ=[];
    end;
else
    J=idenObj(useq1,tseq1,tt,history,state_system_u,pro_system,[]);
end;
function [J,dJ,tt,xx,ttlam,clam,grad]=idenObj(useq,tseq,tt,history,state_system,pro_system,measure_system)
f=state_system.f;fx=state_system.fx;fxd=state_system.fxd;
fu=state_system.fu;fud=state_system.fud;dimx=state_system.dimx;
tao=state_system.tao;
taou=state_system.taou;
if isstruct(pro_system)%如果需要考虑过程目标
    L=pro_system.L;Lx=pro_system.Lx;Lu=pro_system.Lu;
end;
% if isstruct(measure_system);%如果需要考虑辨识目标
% hx=measure_system.hx;
% tsample=measure_system.tsample;
% ysample=measure_system.ysample;
% invQ=measure_system.invQ;
% end;
T=tt(end);
dimtao=numel(tao);dimu=size(useq,1);
tt=unique(union(history.t,tt));
xx=[history.x,zeros(dimx,numel(tt(tt>history.t(end))))];
n=numel(tt);
start=find(tt==history.t(end));
if nargout>1;
    ttlam=[];cfx=zeros(dimx*dimx,n-start+1);cfxd=zeros(dimx*dimx,dimtao,n-start+1);
    cfu=zeros(dimx*dimu,n-start+1);cfud=zeros(dimx*dimu,n-start+1);ulam=zeros(dimu,n-start+1);
end;
ud=[];%如果没有控制时滞，则ud保持这个值
for i=start:n-1;
    x=xx(:,i);
    t=tt(i);h=tt(i+1)-t;
    xh1=interp1(tt,xx',t-tao,'linear','extrap')'; %计算t-lags时刻的状态向量
    u=getu(t,useq,tseq);
    if numel(taou)>0;ud=getu(t-taou,useq,tseq);end
    k1=f(t,x,xh1,u,ud);k1=k1(:);%计算t时刻的dx/dt
    if nargout>1
        ttlam=[ttlam,t];
        fx1=fx(t,x,xh1,u,ud);fxd1=fxd(t,x,xh1,u,ud);fu1=fu(t,x,xh1,u,ud);
        if numel(taou)>0;fud1=fud(t,x,xh1,u,ud);end;
        cfx(:,i-start+1)=fx1(:);
        cfxd(:,:,i-start+1)=reshape(fxd1,dimx*dimx,dimtao);
        cfu(:,i-start+1)=fu1(:);
        if numel(taou)>0;cfud(:,i-start+1)=fud1(:);end;
        ulam(:,i-start+1)=u(:);
    end
    xx(:,i+1)=x(:)+k1(:)*h;%粗预测t+h时刻的状态向量，这样对xx给出t+h时刻的值以后，可以有h>lags
    xh2=interp1(tt,xx',t-tao+h,'linear','extrap')';%插值计算t-lags+h时刻的状态向量，
    k2=f(t,x+h*k1,xh2,u,ud);k2=k2(:);
    x=x+h/2*(k1+k2);
    xx(:,i+1)=x(:);
end;
if nargout>1;
    ttlam=[ttlam,tt(end)];
    fx1=fx(t,x,xh2,u,ud);fxd1=fxd(t,x,xh2,u,ud);fu1=fu(t,x,xh2,u,ud);
    if numel(taou)>0;fud1=fud(t,x,xh2,u,ud);end;
    cfx(:,n-start+1)=fx1(:);cfxd(:,:,n-start+1)=reshape(fxd1,dimx*dimx,dimtao);
    cfu(:,n-start+1)=fu1(:);
    if numel(taou)>0;cfud(:,n-start+1)=fud1(:);end;
    ulam(:,n-start+1)=useq(:,end);
    for i=1:dimx;xxlam(i,:)=interp1(tt,xx(i,:),ttlam,'linear','extrap');end
end;
%  figure;plot(tt,xx);
%过程代价
li=[];
for i=1:numel(tt);
    x=xx(:,i);t=tt(i);
    if t>=0;
        u=getu(t,useq,tseq);
        li=[li,L(t,x,u)];
    end
end
J=trapz(tt(tt>=0),li);
lam0=zeros(dimx,1);%
% J=J+1e6*xx(1,end)+1e6*xx(2,end);lam0=[1e6;1e6];
if nargout>1;
    %将cfxd转成正向时滞
    cfxdnew=0*cfxd;
    for i=1:dimtao;
        li=ttlam+tao(i);
        cfxdi=reshape(cfxd(:,i,:),size(cfxd,1),size(cfxd,3));
        cfxdnew(:,i,li<=T)=interp1(ttlam,cfxdi',li(li<=T))';%,'linear','extrap')';
    end;
    if numel(taou)>0;
        li=ttlam+taou;cfudnew=0*cfud;
        cfudnew(:,li<=T)=interp1(ttlam,cfud',li(li<=T))';%,'linear','extrap')';
    end
    %求解协态方程
    clam=zeros(dimx,numel(ttlam)); lam=lam0;
    grad=zeros(dimu,numel(ttlam));
    %     if isstruct(measure_system);
    %          ksam=numel(tsample);
    %     end;
    for i=1:numel(ttlam)-1;
        t=ttlam(end-i+1);t1=ttlam(end-i);h=t-t1;
        if numel(taou)>0;
            if t+taou<=T;
                lamdu=interp1(ttlam,clam',t+taou)';%,'linear','extrap')';
            else lamdu=zeros(dimx,1);
            end
        end;
        fxi=reshape(cfx(:,end-i+1),dimx,dimx);fxi1=reshape(cfx(:,end-i),dimx,dimx);
        fxdi=reshape(cfxdnew(:,:,end-i+1),dimx,dimx,dimtao);fxdi1=reshape(cfxdnew(:,:,end-i),dimx,dimx,dimtao);
        x=xxlam(:,end-i+1);x1=xxlam(:,end-i);u=ulam(:,end-i+1);u1=ulam(:,end-i);
        fui=reshape(cfu(:,end-i+1),dimx,dimu);
        if numel(taou)>0;fudi=reshape(cfudnew(:,end-i+1),dimx,dimu);else fudi=0;end;
        gradi=fui'*lam+(Lu(t,x,u)+Lu(t1,x1,u1))/2;
        if numel(taou)>0; if t+taou<=T;gradi=gradi+sum(fudi'*lamdu,2);end;end;
        grad(:,end-i+1)=gradi(:);
        Lxi=Lx(t,x,u);Lxi1=Lx(t1,x1,u1);
        if t+tao<=T;
            lamd=interp1(ttlam,clam',t+tao)';%,'linear','extrap')';
        else lamd=zeros(dimx,dimtao);
        end;
        %dlam=Lxi+fxi'*lam+mul(fxdi,lamd,[1,3],[1,2]);
        dlam=Lxi+fxi'*lam+mul_transpose(fxdi,lamd);
        lam1=lam+h*dlam(:);
        clam(:,end-i)=lam1;
        if t1+tao<=T;lamd1=interp1(ttlam,clam',t1+tao)';%,'linear','extrap')';
        else lamd1=zeros(dimx,dimtao);end;
        %dlam1=Lxi1+fxi1'*lam1+mul(fxdi1,lamd1,[1,3],[1,2]);
        dlam1=Lxi1+fxi1'*lam1+mul_transpose(fxdi1,lamd1);
        lam=lam+h/2*(dlam+dlam1);
        clam(:,end-i)=lam(:);
    end;
    x=xxlam(:,1);u=ulam(:,1);
    fui=reshape(cfu(:,1),dimx,dimu);
    if numel(taou)>0;fudi=reshape(cfudnew(:,1),dimx,dimu);else fudi=0;end;
    gradi=fui'*lam+Lu(t,x,u);
    if numel(taou)>0;if t+taou<=T;gradi=gradi+sum(fudi'*lamdu,2);end;end;
    grad(:,1)=gradi;
    
    dJ=zeros(size(useq));
    for i=1:size(useq,2);
        flag=(ttlam<=tseq(i+1)&ttlam>=tseq(i));
        tt_ui=ttlam(flag==1);
        dJ(:,i)=trapz(tt_ui,grad(:,flag==1),2);
    end;
end;

function sumc=mul_transpose(fxd,lamd)
%将fxd中前两个维度先转置，然后与lamd乘
%lamd是dimx,dimtao矩阵；fxd可能是dimx,dimx,dimtao矩阵或者dimx,dimx矩
[dimf,dimx,dimtao]=size(fxd);
if size(lamd,1)~=dimf;
    error;
end
c=zeros(dimx,dimtao);
for k=1:dimtao;
    ak=fxd(:,:,k);
    bk=lamd(:,k);
    c(:,k)=ak'*bk;
end
sumc=sum(c,2);

function u=getu(t,useq,tseq)
if numel(t)==1;
    i=dingwei(t,tseq);
    if i==size(useq,2)+1;u=useq(:,end);
    else u=useq(:,i);
    end;
else
    for k=1:numel(t);
        i=dingwei(t(k),tseq);
           if i==size(useq,2)+1;u=useq(:,end);
            else u=useq(:,i);
            end;
    end;
end;
function i=dingwei(t,tseq);
n=numel(tseq);
if t>tseq(end)|t<tseq(1);'out of range',wrong;
elseif t==tseq(1);i=1;
elseif t==tseq(end);i=n;
else
    i=floor(n/2);ii=[0,n];
    while (t<tseq(i))|t>=tseq(i+1)
        if t<tseq(i);
            ii(2)=i;
            i=floor(mean(ii));
        else
            ii(1)=i+1;
            i=floor(mean(ii));
        end;
        %         i,pause(0.3)
    end;
end
function [J,dJ]=infoOptim(u,tseq,rdfun,parameter,xhistory,uhistory)
useq=reshape(u,2,numel(u)/2);
if nargout>1;
    [J,dJ]=ObjForInfoOptim(useq,tseq,rdfun,parameter,xhistory,uhistory);
else
    J=ObjForInfoOptim(useq,tseq,rdfun,parameter,xhistory,uhistory);
end
K=1e0;J=J*K;
if nargout>1;
    dJ=dJ(:);dJ=dJ*K;dJ=dJ(1:end);
end
function [J,dJ]=infoOptim_x0(u_x,tseq,rdfun,parameter,uhistory)
%考虑了初始值的优化
global dimx;
xhistory=u_x(1:dimx);
u=u_x(dimx+1:end);
useq=reshape(u,2,numel(u)/2);
if nargout>1;
    [J,dJ,dJx0]=ObjForInfoOptim(useq,tseq,rdfun,parameter,xhistory,uhistory);
else
    J=ObjForInfoOptim(useq,tseq,rdfun,parameter,xhistory,uhistory);
end
K=1e0;J=J*K;
if nargout>1; dJ=[dJx0(:);dJ(:)];dJ=dJ*K;end