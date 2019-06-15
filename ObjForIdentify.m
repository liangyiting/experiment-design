function [J,dJp,tt,zz]=ObjForIdentify(rdseq,tseq_rd,parameter,useq,tseq,xhistory,uhistory,measure_system)
global dimtheta dimtao dimx 
wan_ge_mi_du=40;%为计算信息判据积分系统时所采用的网格密度%计算网格不能小于20.Ngrid_越大，梯度计算越准确 
if nargout>-1;
    rdfun=@(t)interp_piecewise(tseq_rd,rdseq,t);
    [zz,tt]=ObjForInfoOptim(useq,tseq,rdfun,parameter,wan_ge_mi_du,xhistory,uhistory,measure_system);
    tsample=measure_system.tsample;ysample=measure_system.ysample;hx=measure_system.hx;invQ=measure_system.invQ;
    ysimulation=hx*interp1(tt,zz(1:dimx,:)',tsample)'; 
    sen=interp1(tt,zz(dimx+1:end,:)',tsample)';
    nsample=numel(tsample);
    dxtheta=reshape(sen(1:dimx*dimtheta,:),dimx,dimtheta,nsample);
    dxtao=reshape(sen(dimx*dimtheta+1:dimx*dimtheta+dimx*dimtao,:),dimx,dimtao,nsample);
    dJtao=zeros(1,dimtao);dJtheta=zeros(1,dimtheta);
    J=0;cj=[];
    for i=1:nsample;
        er=ysimulation(:,i)-ysample(:,i);
        dJtao=dJtao+1*er(:)'*invQ*hx*dxtao(:,:,i);
        dJtheta=dJtheta+1*er(:)'*invQ*hx*dxtheta(:,:,i);
        J=J+1/2*er(:)'*invQ*er(:);cj=[cj,J];
    end;
    dJp=[dJtheta(:);dJtao(:)]';
end;
function [zz,tt]=ObjForInfoOptim(useq,tseq,rdfun,parameter,wan_ge_mi_du,xhistory,uhistory,measure_system)
% global augmented_system imp_system measure_system;
global g gz gzd gu gud taoz taouz
global f fx fxd fu fud fr tao taou theta;
global impfun impfz impfzd impfu impfud taoimp taouimp;
global dimtheta dimtao dimx;
xinitial.x=[xhistory(:),xhistory(:)];
xinitial.t=[-50,0];%足够负，至少要比tao允许的最大负值要负
theta=parameter(1:dimtheta);tao=parameter(dimtheta+1:dimtheta+dimtao);taou=parameter(dimtheta+dimtao+1:end);
uinitial.u=uhistory(:);
dimz=dimx+dimx*(dimtheta+dimtao+1);
history.x=[xinitial.x;zeros(dimz-dimx,numel(xinitial.t))];
history.t=xinitial.t;
history.u=uinitial.u;
% history.xt=0*history.x;
T=tseq(end);
t_interval=[0,T];
chix=@(t)(t-tao>=0);
augmented_system=struct('f',@(t,x,xd,u,ud)g(x,xd,u,ud,chix(t),rdfun(t),parameter),...
    'fx',@(t,x,xd,u,ud)gz(x,xd,u,ud,chix(t),rdfun(t),parameter),...
    'fxd',@(t,x,xd,u,ud)gzd(x,xd,u,ud,chix(t),rdfun(t),parameter),...
    'fu',@(t,x,xd,u,ud)gu(x,xd,u,ud,chix(t),rdfun(t),parameter),...
    'fud',@(t,x,xd,u,ud)gud(x,xd,u,ud,chix(t),rdfun(t),parameter),...
    'tao',taoz(parameter),'taou',taouz(parameter),'dimx',dimz);
parameter=struct('tao',tao,'taou',taou,'theta',theta);
[zz,tt]=Obj(useq,tseq,t_interval,history,augmented_system,measure_system,wan_ge_mi_du);
function [zz,tt]=Obj(useq,tseq,t_interval,history,state_system,measure_system,wan_ge_mi_du)
%1. 求解灵敏度方程 ---以及伴随系统--以及梯度
g=state_system.f;
taoz=state_system.tao;
taouz=state_system.taou;
% hx=measure_system.hx;
tsample=measure_system.tsample;
% dimu=size(useq,1);
% dimtaoz=numel(taoz);dimtaouz=numel(taouz);
%在<0时刻用历史的u来计算
useq=[history.u,useq];
tseq=[history.t(1:end-1),tseq];

%给出合适的计算网格
Ngrid_=floor(wan_ge_mi_du*(t_interval(end)-t_interval(1)));
tt=genGrid(tsample,tseq,t_interval,Ngrid_);
history=update(history,g,tt,tseq,useq,taoz,taouz);
zz=history.x;tt=history.t;
function [tt,tpulse]=genGrid(tsample,tseq,t_interval,Ngrid)
tstart=t_interval(1);T=t_interval(2);
tt=linspace(tstart,T,Ngrid);
tpulse=union(tseq,tsample);ep=1e-4;
tt=sort(union(tt,[tpulse,tpulse-ep,tpulse+ep]));
tt=tt(tt>=tstart&tt<=T);
function history= update(history,g,tt,tseq,useq,tao,taou)
ddefun=@(t,x,xh,u,ud)g(t,x,xh,u,ud);
x=history.x(:,end);start=numel(history.t);
N=numel(tt);
for i=1:N;
    t=tt(i);
    xh1=interp1(history.t,history.x',t-tao)';%计算t-lags时刻的状态向量
    u=getu(t,useq,tseq);ud=getu(t-taou,useq,tseq);
    k1=ddefun(t,x,xh1,u,ud);k1=k1(:);%计算t时刻的dx/dt
    if i+1<=N;
        h=tt(i+1)-t;
        history.t(i+start)=tt(i+1);    
        history.x(:,i+start)=x(:)+k1(:)*h;%粗预测t+h时刻的状态向量，这样对xx给出t+h时刻的值以后，可以有h>lags 
        xh2=interp1(history.t,history.x',t-tao+h)';%插值计算t-lags+h时刻的状态向量，
        try
        u=getu(tt(i+1),useq,tseq);
        catch
            1;
        end
       ud=getu(tt(i+1)-taou,useq,tseq); 
        k2=ddefun(t+h,x+h*k1,xh2,u,ud);k2=k2(:);
        x=x+h/2*(k1+k2);
        history.x(:,i+start)=x; 
    end
end;
function u=getu(t,useq,tseq)
if numel(t)==1;
    i=dingwei(t,tseq);
    if i==size(useq,2)+1;u=useq(:,end);
    else u=useq(:,i);
    end;
else
    u=[]; for k=1:numel(t);u(:,k)=getu(t(k),useq,tseq);end;
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