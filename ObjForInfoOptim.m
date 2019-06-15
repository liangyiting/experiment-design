function [J,dJ,dJx0,newhistory,H]=ObjForInfoOptim(useq,tseq,rdfun,parameter,xhistory,uhistory)
% global augmented_system imp_system measure_system;
global g gz gzd gu gud taoz taouz
global f fx fxd fu fud fr tao taou theta;
global impfun impfz impfzd impfu impfud taoimp taouimp;
global measure_system;
global dimtheta dimtao dimx dimu;
xinitial.x=[xhistory(:),xhistory(:)];
xinitial.t=[-20,0];
% dimtheta=2;dimtao=2;dimx=2;dimu=2;
dimz=dimx*(1+numel(parameter));
theta=parameter(1:dimtheta);tao=parameter(dimtheta+1:dimtheta+dimtao);taou=parameter(dimtheta+dimtao+1:end);
xthetaInitial=zeros(1,dimtheta);xtaoInitial=zeros(1,dimtao);xtaouInitial=0;
uinitial.u=uhistory(:);
dimx=2;dimz=dimx+dimx*(dimtheta+dimtao+1);
history.x=[xinitial.x;zeros(dimz-dimx,numel(xinitial.t))];
history.t=xinitial.t;
history.u=uinitial.u;
history.xt=0*history.x;
T=tseq(end);
t_interval=[0,T];
chix=@(t)(t-tao>=0);
augmented_system=struct('f',@(t,x,xd,u,ud)g(x,xd,u,ud,chix(t),rdfun(t),parameter),...
    'fx',@(t,x,xd,u,ud)gz(x,xd,u,ud,chix(t),rdfun(t),parameter),...
    'fxd',@(t,x,xd,u,ud)gzd(x,xd,u,ud,chix(t),rdfun(t),parameter),...
    'fu',@(t,x,xd,u,ud)gu(x,xd,u,ud,chix(t),rdfun(t),parameter),...
    'fud',@(t,x,xd,u,ud)gud(x,xd,u,ud,chix(t),rdfun(t),parameter),...
    'tao',taoz(parameter),'taou',taouz(parameter),'dimx',dimz);
imp_system=struct('f',@(t,x,xd,u,ud)impfun(x,xd,u,ud,rdfun(t),parameter),...
    'fx',@(t,x,xd,u,ud)impfz(x,xd,u,ud,rdfun(t),parameter),...
    'fxd',@(t,x,xd,u,ud)impfzd(x,xd,u,ud,rdfun(t),parameter),...
    'fu',@(t,x,xd,u,ud)impfu(x,xd,u,ud,rdfun(t),parameter),...
    'fud',@(t,x,xd,u,ud)impfud(x,xd,u,ud,rdfun(t),parameter),...
    'tao',taoimp(parameter),'taou',taouimp(parameter));
parameter=struct('tao',tao,'taou',taou,'theta',theta);
wan_ge_mi_du=40;%����ģ�������õ������ܶ�
if nargout>1;
    [J,dJ,dJx0,newhistory,~,H]=Obj(useq,tseq,parameter,t_interval,history,augmented_system,imp_system,measure_system,wan_ge_mi_du);
else
    J=Obj(useq,tseq,parameter,t_interval,history,augmented_system,imp_system,measure_system,wan_ge_mi_du);
end;
function [J,dJ,dJx0,history,tt,H,Hs]=Obj(useq,tseq,parameter,t_interval,history,state_system,imp_system,measure_system,wan_ge_mi_du)
T=t_interval(2);tstart=t_interval(1);%�������
Ngrid_=(T-tstart)*wan_ge_mi_du;%�Զ��������ܶȣ���Ҫ
try tao=parameter.tao;dimtao=numel(tao);end;
try taou=parameter.taou;end;
try theta=parameter.theta;dimtheta=numel(theta);end;
%��Ҫ������������history.x��history.t,�Լ�history.u������˵��ֱ�������е�history
%1. ��������ȷ��� ---�Լ�����ϵͳ--�Լ��ݶ�
global zlaminitial ;
g=state_system.f;gz=state_system.fx;gzd=state_system.fxd;
gu=state_system.fu;gud=state_system.fud;
dimz=state_system.dimx;
taoz=state_system.tao;
taouz=state_system.taou;

impfun=imp_system.f;
impfz=imp_system.fx;
impfzd=imp_system.fxd;
impfu=imp_system.fu;
impfud=imp_system.fud;
taoimp=imp_system.tao;
taouimp=imp_system.taou;

hx=measure_system.hx;
tsample=measure_system.tsample;
invQ=measure_system.invQ;

dimu=size(useq,1);Nu=size(useq,2);
dimtaoz=numel(taoz);dimtaouz=numel(taouz);
dimx=size(hx,2);

%��<0ʱ������ʷ��u������
numel_historyu=size(history.u,2);
useq=[history.u,useq];
tseq=[history.t(1:end-1),tseq];

%�������ʵļ�������
[tt,tpulse,tpulselam,tpulse_1,tpulse_2,tpulse_3]=genGrid(tsample,tseq,t_interval,tao,taou,taoz,taouz,taoimp,taouimp,Ngrid_);
%������ʷ�����㷽���ݶ���Ϣ
info.fun={gz,gzd,gu,gud};
info.funeval={zeros(dimz*dimz,2),zeros(dimz^2*dimtaoz,2),zeros(dimz*dimu,2),zeros(dimz*dimtaouz*dimu,2)};
%
infoimp=struct('fz',{{}},'fzd',{{}},'fu',{{}},'fud',{{}},'zd',{{}},'ud',{{}});
ttk=tt(tt<=tpulse(1)&tt>=tstart);
%%
%1. ���
[history,info]=update(history,g,ttk,tseq,useq,taoz,taouz,info);
timp=zeros(1,numel(tpulse));
for k=1:numel(tpulse);
    t=tpulse(k);if k~=numel(tpulse);tnext=tpulse(k+1);else tnext=T;end;
    ep=1e-7;uall=getu([t,t-taouimp(:)'],useq,tseq);u=uall(:,1);udimp=uall(:,2:end);
    zdimp=interp1(history.t,history.x',t-taoimp)';
    infoimp.zd{k}=zdimp;infoimp.ud{k}=udimp;
    z=history.x(:,end);
    znew=impfun(t,z,zdimp,u,udimp);
    
    timp(k)=t;
    if ismember(t,tseq+taou);
        infoimp.fz{k}=reshape(impfz(t,z,zdimp,u,udimp),dimz,dimz);
        infoimp.fzd{k}=reshape(impfzd(t,z,zdimp,u,udimp),dimz,dimz,numel(taoimp));
        infoimp.fu{k}=reshape(impfu(t,z,zdimp,u,udimp),dimz,dimu);
        infoimp.fud{k}=reshape(impfud(t,z,zdimp,u,udimp),dimz,dimu,numel(taouimp));
    end;
    m=numel(history.t);
    %     history.t(m:m+1)=[t-ep,t];
    history.t(m:m+1)=[1/2*(history.t(m)+history.t(m-1)),t];
    history.x(:,m+1)=znew;
    ttk=tt(tt>=t&tt<=tnext);[history,info]=update(history,g,ttk,tseq,useq,taoz,taouz,info);
end;
zz=history.x;tt=history.t;
% figure;plot(tt,zz(1:2,:)');
%%
Lz=zeros(dimz,numel(tt));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dimtheta=numel(theta);
dimp=dimtao+dimtheta+1;
ztsam=interp1(tt,zz',tsample)';
H=zeros(dimp);
cM={};
for i=1:length(tsample);
    xp=reshape(ztsam(dimx+1:end,i),dimx,dimp);
    M=xp'*hx'*invQ*hx;H=H+M*xp;
    cM{i}=M;
end
%1. ��Ϣָ��Ͷ�Ӧ����Ϣָ��J�Ը���ʱ�̵�״̬�������ݶ�dJz
global criterion;
if numel(criterion)==0; criterion='D';end
dJz=zeros(dimz,numel(tsample));
global len_p;%
s=len_p;
%s=dimtheta+dimtao;
for i=1:length(tsample);
    M=cM{i};
    switch criterion
        case 'D'
            %C=inv(H);J=-det(H);dJz(dimx+1:end,i)=-2*kron(eye(dimp),M)'*(det(H)*C(:));
            C=inv(H(1:s,1:s));J=-det(H(1:s,1:s));dJz(dimx+1:dimx+s*dimx,i)=-2*kron(eye(s),M(1:s,:))'*(-J*C(:));
        case 'A'
            %C=inv(H);J=trace(C);ctc=C'*C;dJz(dimx+1:end,i)=-2*kron(eye(dimp),M)'*ctc(:);
            Hs=H(1:s,1:s);C=inv(Hs);J=trace(C);ctc=C'*C;dJz(dimx+1:dimx+dimx*s,i)=-2*kron(eye(s),M(1:s,:))'*ctc(:);
        case 'E'
            Hs=H(1:s,1:s); [v,D]=eig(Hs);
            d=diag(D);[dnew,index]=sort(d);imin=1;J=-dnew(imin); vmin=v(:,index(imin));V=vmin*vmin';
            dJz(dimx+1:dimx+dimx*s,i)=-2*kron(eye(s),M(1:s,:))'*V(:);
        case 'vA'
            %J=-trace(H);ctc=eye(dimp);dJz(dimx+1:end,i)=-2*kron(eye(dimp),M)'*ctc(:);
            J=-trace(H(1:s,1:s));ctc=eye(s);dJz(dimx+1:dimx+dimx*s,i)=-2*kron(eye(s),M(1:s,:))'*ctc(:);
    end;
end
%2. ���̾�����ָ��Ͷ�״̬���ݶ�
%3. ����Լ���Ͷ�״̬���ݶ�
%2.��ʶָ��ͱ�ʶָ���״̬���ݶ�
test=0;if test==1;Kq=2;J=sum(ztsam(Kq,:));end;
%%
if nargout==1;return;end;%�����Ҫ�����ݶ�
%ת��Ϊlam��������ʱ
cgz=reshape(info.funeval{1},dimz,dimz,numel(tt));
cgu=reshape(info.funeval{3},dimz,dimu,numel(tt));
cgzd=reshape(info.funeval{2},dimz*dimz,dimtaoz,numel(tt));
cgud=reshape(info.funeval{4},dimz*dimu,dimtaouz,numel(tt));
cgzdnew=zeros(dimz,dimz,dimtaoz,numel(tt));cgudnew=zeros(dimz,dimu,dimtaouz,numel(tt));
for j=1:dimtaoz;
    ttj=tt+taoz(j);
    cgzdnew(:,:,j,ttj<=T)=reshape(interp1(tt,reshape(cgzd(:,j,:),dimz*dimz,numel(tt))',ttj(ttj<=T))',dimz,dimz,sum(ttj<=T));
end;
for j=1:dimtaouz;
    ttj=tt+taouz(j);
    cgudnew(:,:,j,ttj<=T)=reshape(interp1(tt,reshape(cgud(:,j,:),dimz*dimu,numel(tt))',ttj(ttj<=T))',dimz,dimu,sum(ttj<=T));
end;
%%%%%%%%%%%%%%%%%%%%%%%%
lam=zlaminitial; n=numel(tt);
history.lam(:,n)=lam(:);
grad=zeros(dimu,n);

isjump{1}=ismember(tt,tpulse_1);%��һ������
isjump{2}=ismember(tt,tpulse_2);%�ڶ�������
for i=1:numel(tpulse_3);isjump{i+2}=ismember(tt,tpulse_3{i});end;%����������
grad=zeros(dimu,n);
for i=1:n;
    t=history.t(n-i+1);lam=history.lam(:,n-i+1);lam0=lam;
    %         %������巽��
    if isjump{1}(n-i+1);
        j=find(tsample==t);
        if test==0;
            lam=lam+dJz(:,j);
        else
            li=zeros(dimz,1);li(Kq)=1;lam=lam+li;%%�������ȵĵ�Kq�����������Լ���׼ȷ��
        end;
    end;%��Ӧ��
    
    
    if isjump{2}(n-i+1);%
        k=find(timp==t);
        %%
        lam=infoimp.fz{k}'*lam0-lam0+lam;%+lam;
    end
    %%%%%%%%%%%%%%%%
    %
    for j=1:numel(tpulse_3);
        if isjump{j+2}(n-i+1);
            [~,k]=min(abs(timp-t-tao(j)));
            [~,s]=min(abs(tt-t-tao(j)));
            lamdimpj=history.lam(:,s);
            %%
            lam=infoimp.fzd{k}(:,:,j)'*lamdimpj+lam;
        end
    end;
    
    
    %����Ŀ�꺯�����ݶ�
    gui=cgu(:,:,n-i+1);
    gudi=cgudnew(:,:,:,n-i+1);
    lamdu=zeros(dimz,numel(taouz));li=t+taouz;
    lamdu(:,li<=T)=interp1(tt,history.lam',li(li<=T),'nearest')';
    grad(:,n-i+1)=grad(:,n-i+1)+gui'*lam+mul(gudi,lamdu,[1,3],[1,2]);%%���Ĺ�ʽ
    %����lam
    if (i<n);
        %��ȡ�����״̬���̵�����Ϣ
        gzi=cgz(:,:,n-i+1);
        gzdi=cgzdnew(:,:,:,n-i+1);
        lamd=zeros(dimz,dimtaoz);s=t+taoz(:);lamd(:,s<=T)=interp1(tt,history.lam',s(s<=T),'nearest')';
        k1=Lz(:,n-i+1)+gzi'*lam+mul(gzdi,lamd,[1,3],[1,2]);
        h=t-history.t(n-i);
        lamnext=lam+k1(:)*h;history.lam(:,n-i)=lamnext;
        s=t+taoz(:);lamd(:,s<=T)=interp1(tt,history.lam',s(s<=T)-h,'nearest')';
        gzi=cgz(:,:,n-i);gzdi=cgzdnew(:,:,:,n-i);
        k2=Lz(:,n-i)+gzi'*lamnext+mul(gzdi,lamd,[1,3],[1,2]);%%���Ĺ�ʽ
        lam=lam+h/2*(k1+k2);
        history.lam(:,n-i)=lam;
    end;
end;
dJx0=lam;
%�ۺϵó���Ŀ������ĵ���
dJ=zeros(dimu,Nu);
for i=numel_historyu+1:size(useq,2);%useq�ĵ�һ������ʷ���룬�ڶ��ʼ������Ƶ�����
    flag=(tt<=tseq(i+1)&tt>=tseq(i));
    tt_ui=tt(flag==1);
    dJ(:,i-1)=dJ(:,i-1)+trapz(tt_ui,grad(:,flag==1),2);
end;
%%
%�����ȷ��̵����巽����uӰ����z�����Ӱ����J�������֣�u(t-ep),u(t+ep)
for k=1:numel(tpulse_2);
    %%grad�ļ�����û�п���infoimp.fud��infoimp.fu!��
    %%%%%
    t=tpulse_2(k);
    lamk=interp1(tt,history.lam',t-ep,'nearest','extrap')';%g/z(t)
    %[~,k1]=min(abs(tseq(numel_historyu+1:end)-(t-taou)));
    k1=k;
    dJ(:,k1)=dJ(:,k1)+infoimp.fud{k}(:,:,1)'*lamk;%g/u(t-(taou-ep))����unew�ĵ���
    if k1>1;
        dJ(:,k1-1)=dJ(:,k1-1)+infoimp.fud{k}(:,:,2)'*lamk;%g/u(t-(taou+ep)),��uold�ĵ���
    end;
end;

if nargout>2;
    tt1=linspace(0,max(taoz),20);
    grad=zeros(dimz,numel(tt1));
    for i=1:numel(tt1);
        t=tt1(i);
        z=interp1(history.t,history.x',t)';
        zd=interp1(history.t,history.x',t-taoz)';
        u=getu(t,useq,tseq);ud=getu(t-taouz,useq,tseq);
        gzdi=reshape(gzd(t,z,zd,u,ud),dimz,dimz,dimtaoz);
        lam=interp1(tt,history.lam',t)';lam=lam(:)*ones(1,numel(taoz));lam(:,t-taoz>0)=0;
        %         lamd=zeros(dimz,dimtaoz);s=t+taoz(:);lamd(:,s<=T)=interp1(tt,history.lam',s(s<=T),'nearest')'
        grad(:,i)=mul(gzdi,lam,[1,3],[1,2]);
    end
    dJx0=interp1(tt,history.lam(1:2,:)',0)'+trapz(tt1,grad(1:dimx,:),2);
end
function [tt,tpulse,tpulselam,tpulse_1,tpulse_2,tpulse_3]=genGrid(tsample,tseq,t_interval,tao,taou,taoz,taouz,taoimp,taouimp,Ngrid)
tstart=t_interval(1);
T=t_interval(2);
dimtao=numel(tao);
tpulse_1=tsample;%Э̬���̵�һ�������-��ΪĿ�꺯���д��ڼ�ϵ�
tpulse_2=taou+tseq;tpulse_2=tpulse_2(tpulse_2>tstart&tpulse_2<T);%Э̬���̵ڶ��������--��Ϊģ���д��ڿ���ʱ��
tpulse=[tpulse_2];%�����ȷ��������
tpulse_3={};%for i=1:dimtao;li=tpulse_2-tao(i);tpulse_3{i}=li(li>=tstart);end;%Э̬���̵����������
tpulselam=[tpulse_1,tpulse_2];for i=1:numel(tpulse_3);tpulselam=[tpulselam,tpulse_3{i}];end;%���е������
tpulseTaouimp=[];
ep=1e-4;
try Ngrid;catch; Ngrid=40*(T-tstart);end;
%tt=linspace(tstart-5,T,Ngrid);
tt=unique([linspace(tstart-5,0,10),linspace(0,T,Ngrid)]);
tt=sort(union(tt,[tpulselam,tpulselam-ep,tpulselam+ep,tpulseTaouimp]));
function [history,info]= update(history,g,tt,tseq,useq,tao,taou,info,rdfun)
try rdfun(tt(1));
    ddefun=@(t,x,xh,u,ud)g(t,x,xh,u,ud,rdfun(t));
catch rdflag=0;
    ddefun=@(t,x,xh,u,ud)g(t,x,xh,u,ud);
end;
x=history.x(:,end);start=numel(history.t);
N=numel(tt);
for i=1:N;
    t=tt(i);
    xh1=interp1(history.t,history.x',t-tao)';%����t-lagsʱ�̵�״̬����
    u=getu(t,useq,tseq);ud=getu(t-taou,useq,tseq);
    if nargout>1;
        for j=1:numel(info.fun);
            li=info.fun{j}(t,x,xh1,u,ud);
            info.funeval{j}(:,i+start-1)=li(:);
        end;
    end;
    k1=ddefun(t,x,xh1,u,ud);k1=k1(:);%����tʱ�̵�dx/dt
    if i+1<=N;
        h=tt(i+1)-t;
        history.t(i+start)=tt(i+1);    
        history.x(:,i+start)=x(:)+k1(:)*h;%��Ԥ��t+hʱ�̵�״̬������������xx����t+hʱ�̵�ֵ�Ժ󣬿�����h>lags 
        xh2=interp1(history.t,history.x',t-tao+h)';%��ֵ����t-lags+hʱ�̵�״̬������
        u=getu(tt(i+1),useq,tseq);
       try ud=getu(tt(i+1)-taou,useq,tseq); catch
           1
       end
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
    for k=1:numel(t);
        %i=dingwei(t(k),tseq);
        u(:,k)=getu(t(k),useq,tseq);
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


function sumc=mul(fxd,lamd,ind1,ind2)
%function sumc=mul_transpose(fxd,lamd)
%����ԭ����mul�е�һ������ԭ����mul��dimtao>1ʱ�ǶԵģ�����=1�ǻ��������Ľ��
%��fxd��ǰ����ά����ת�ã�Ȼ����lamd��
%lamd��dimx,dimtao����fxd������dimx,dimx,dimtao�������dimx,dimx��
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


