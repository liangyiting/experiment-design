function []=getequations()
%һ.����ϵͳģ��
clear all;
model_id=2;
switch model_id 
    case 1
        %ģ��1
        Dim.x=2;Dim.theta=4;Dim.tao=1;Dim.u=2;
        dimx=Dim.x;dimtheta=Dim.theta;dimtao=Dim.tao;dimu=Dim.u; dimrd=2;
        thetaScale=[0.310,0.180,0.550,0.050];
        f=@(x,xd,u,ud,theta,rd)[(theta(1)*thetaScale(1)*x(2)/(theta(2)*thetaScale(2)+x(2))-u(1))*x(1)-theta(4)*thetaScale(4)*xd(1,1);...
            -theta(1)*thetaScale(1)*x(2)/(theta(2)*thetaScale(2)+x(2))*x(1)/theta(3)*thetaScale(3)+u(1)*(u(2)-x(2))];
        savemat_name='model_1.mat';
    case 2 
        %ģ��2
        Dim.x=2;Dim.theta=2;Dim.tao=1;Dim.u=2;
        dimx=Dim.x;dimtheta=Dim.theta;dimtao=Dim.tao;dimu=Dim.u; dimrd=2;
        Q=200;beta=[17,710];V=400;thetaScale=[7.828e-4,2.823e-4];
        f=@(x,xd,u,ud,theta,rd)1/V*[beta(1)*xd(2,1)-Q*xd(1,1)-theta(1)*thetaScale(1)*u(1)*xd(1,1)+Q*rd(1);
            beta(2)*xd(1,1)-Q*xd(2,1)-theta(2)*thetaScale(2)*u(2)*xd(2,1)+Q*rd(2)];
%         f=@(x,xd,u,ud,theta,rd)1/V*[beta(1)*xd(2,1)-Q*xd(1,1)-theta(1)*thetaScale(1)*u(1)*xd(1,1)+Q*rd(1);
%             beta(2)*xd(1,1)-Q*xd(2,1)-theta(2)*1*u(2)*1+Q*rd(2)];
        savemat_name='model_2.mat';
end

%��.���ű�������
%1. ԭϵͳ 
    %1.1����ԭϵͳ����ʶ����
    p=matrix1('p',dimtheta+dimtao+1);
    theta=p(1:dimtheta);tao=p(dimtheta+1:dimtheta+dimtao);taou=p(dimtheta+dimtao+1:end);
    %1.2 ����ԭϵͳ��֪����  
%     syms Q V beta1 beta2 real;% known parameters
    %1.4 ����ʱ�����t
    syms t real;
    %1.5 ����ԭϵͳ״̬�ĳ�ʼֵ-ֻ���ǳ�ʼ״̬Ϊ����������
    xinitial=matrix1('xinitial',dimx);  

%״̬����
    %1 ����ԭϵͳ��״̬�����Ϳ��Ʊ������Ŷ�����
    x=matrix1('x',dimx);
    xd=matrix1('xd',dimx,dimtao);
    u=matrix1('u',dimu);
    ud=matrix1('ud',dimu);
    rd=matrix1('rd',2);%�ⲿ�Ŷ��������ģ���о�����������
    %2 �Ƶ�ԭϵͳЭ̬���������ݶ�
    f_=f(x,xd,u,ud,theta,rd);
    fx_=myjacobian(f_,x);
    fxd_=myjacobian(f_,xd);
    fu_=myjacobian(f_,u);
    fud_=myjacobian(f_,ud);
    fr_=myjacobian(f_,rd);
    
%����ϵͳ    
    %1 ��������ϵͳ��״̬�����ͳ���״̬����
    dimz=dimx*(2+dimtao+dimtheta);
    z=matrix1('z',dimz);
    zd=matrix1('zd',dimz,dimtao+dimtao*dimtao);    
    %2 ��������ϵͳ�Ŀ�������
    u=matrix1('u',dimu);
    uh=matrix1('uh',dimu,1+dimtao*2);
    ud=uh(:,1);utao=uh(:,2:1+dimtao);udtao=uh(:,2+dimtao:2*dimtao+1);    
    %1.6 ����ԭϵͳ��Ҳ������ϵͳ�ģ����ⲿ�Ŷ���δ֪��
    rd=matrix1('rd',dimrd);%�ⲿ�Ŷ��������ģ���о�����������
    %3 ��������ϵͳ��ʱ�Ͳ���
    taouz_=[taou;tao(:);taou+tao(:)];%ud,utao,udtao
    li=[];for i=1:dimtao;li=[li,tao(i)+tao];end;
    taoz_=[tao(:);li(:)];%%state-delay for sensitivity system
    dimtaoz=numel(taoz_);
    dimtaouz=numel(taouz_);
    [x,xtheta,xtao,xtaou]=z_to_x(z,Dim);
    [xd,xthetad,xtaod,xtaoud,xdd]=zd_to_x(zd,Dim);
    
%3. ����ϵͳ
    %4.1 ��������ϵͳ��ʱ�Ͳ���
    taoimp_=tao;
    ep=1e-7;
    taouimp_=[taou-ep,taou+ep];%taou-ep��Ӧ����unew����������ұߵ�����u
    %4.2 ��������ϵͳ�Ŀ��Ʊ���
    udimp=matrix1('udimp',dimu,numel(taouimp_));
    unew=udimp(:,1);
    uold=udimp(:,2);
    %4.3��������ϵͳ��״̬����
    zimp=matrix1('zimp',dimz);
    zdimp=matrix1('zdimp',dimz,dimtao);
    ximp=zimp(1:dimx,:);xdimp=zdimp(1:dimx,:);

%��. �Ƶ�����
%1.����ϵͳģ��
    fv=f(x,xd,u,ud,theta,rd);
    fx=myjacobian(fv,x);
    fxd=myjacobian(fv,xd);
    ftheta=myjacobian(fv,theta);
    %calculate dx/dt(x-tao)=f(t-tao)
    fsym=@(x,i)sym(strcat(x,num2str(i)),'real');
    chi=matrix1('chix',dimtao); %���㷽��ת��ָʾ��
    %������ֻ����xinitialΪ���������
    clear a;clear sdxd;for i=1:dimx;for j=1:dimtao;a(i,j)=sym('0');end;end;sdxd=a;
    for i=1:dimtao;sdxd(:,i)=chi(i).*f(xd(:,i),xdd(:,:,i),utao(:,i),udtao(:,i),theta,rd)+(1-chi(i))*diff(xinitial,t);end;
    for i=1:dimtao;a(:,i)=fxd(:,:,i)*sdxd(:,i);end;
    if dimtao>1;
    dxtheta=fx*xtheta+ftheta+mul(fxd,xthetad,[2,3],[1,3]);
    dxtao=fx*xtao+mul(fxd,xtaod,[2,3],[1,3])-a;
    dxtaou=fx*xtaou+mul(fxd,xtaoud,[2,3],[1,2]);
    elseif dimtao==1;
        dxtheta=fx*xtheta+ftheta+fxd*xthetad;
        dxtao=fx*xtao+fxd*xtaod-a;
        dxtaou=fx*xtaou+fxd*xtaoud;
    end;
    dz=[fv(:);dxtheta(:);dxtao(:);dxtaou(:)];

    %2 ����ϵͳ�ĳ�ʼ״̬
    xthetaInitial=myjacobian(xinitial,theta);
    xtaoInitial=myjacobian(xinitial,tao);
    xtaouInitial=diff(xinitial,taou);
    zinitial=[xinitial(:);xthetaInitial(:);xtaoInitial(:);xtaouInitial(:)]';
%2. �Ƶ�����ϵͳ�ķ���
    xtaouJump=f(ximp,xdimp,u,uold,theta,rd)-f(ximp,xdimp,u,unew,theta,rd);
    zjump=zimp+[zeros(dimx+dimx*dimtheta+dimx*dimtao,1);xtaouJump(:)];
%3. ���㷽��Э̬ϵͳ
    g_=dz;
    %5.2 �Ƶ�Э̬ϵͳ���̵�������������ĸ��ָ��涼
    gz_=myjacobian(g_,z);
    gzd_=myjacobian(g_,zd);
    gu_=myjacobian(g_,u);
    gud_=myjacobian(g_,uh);
    %5.-1 ���巽��ת��ʱ�����
    chilam=matrix1('chilam',dimz,numel(taoz_));
    %5.0 ����Э̬ϵͳ��״̬ʱ�Ͳ���
    taolam_=tao;
    %5.1 ����Э̬����
    lam=matrix1('lam',dimx);
    lamd=matrix1('lamd',dimx,numel(taoz_));
    % 5.3 Э̬�����ĳ�ʼֵ������ֻ����phiΪ�������
    phi=0;
    zlaminitial_=-myjacobian(phi,z);zlaminitial_=zlaminitial_(:);
%6. �Ƶ�Э̬���巽��
    impfz_=myjacobian(zjump,zimp);
    impfzd_=myjacobian(zjump,zdimp);
    impfu_=myjacobian(zjump,u);
    impfud_=myjacobian(zjump,udimp);

% save model_1_symbol;

eval(['taolam=@(p)',changeExp1(taolam_),';']);
%1. ����ϵͳ �� ����ϵͳ ��ʱ�ͣ� �Լ�����ϵͳЭ̬���̵�ĩ��ֵ
eval(['taoz=@(p)',changeExp1(taoz_),';']);
eval(['taouz=@(p)',changeExp1(taouz_),';']); 
eval(['taoimp=@(p)',changeExp1(taoimp_),';']);
eval(['taouimp=@(p)',changeExp1(taouimp_),';']);
eval(['zlaminitial=',changeExp1(zlaminitial_),';']);
%2. ����ϵͳ������g�͵���gz,gzd,...
% sg=g;sgz=gz;sgzd=gzd;sgu=gu;sgud=gud;%������ű��ʽ
eval(['g=@(z,zd,u,uh,chix,rd,p)',changeExp1(g_),';']);
eval(['gz=@(z,zd,u,uh,chix,rd,p)',changeExp1(gz_),';']);
eval(['gzd=@(z,zd,u,uh,chix,rd,p)',changeExp1(gzd_),';']);
eval(['gu=@(z,zd,u,uh,chix,rd,p)',changeExp1(gu_),';']);
eval(['gud=@(z,zd,u,uh,chix,rd,p)',changeExp1(gud_),';']);
% 3. ���巽�̺��䵼�� impfun,impfz,impfzd,impfu,impfud.
eval(['impfun=@(zimp,zdimp,u,udimp,rd,p)',changeExp1(zjump),';']);%���巽��
eval(['impfz=@(zimp,zdimp,u,udimp,rd,p)',changeExp1(impfz_),';']);%���巽����
eval(['impfzd=@(zimp,zdimp,u,udimp,rd,p)',changeExp1(impfzd_),';']);
eval(['impfu=@(zimp,zdimp,u,udimp,rd,p)',changeExp1(impfu_),';']);
eval(['impfud=@(zimp,zdimp,u,udimp,rd,p)',changeExp1(impfud_),';']);
%2.1 ԭϵͳ����f,fx,fxd,fu,fud,fr
% gen_stateEq_model_1;
% parameters_case_model_1;%����ģ���а����˲��ַǱ�ʶ������Q,V�ȣ�ÿһ�η������㶼�����°���Щ��������ֵΪ������������֮��Ҫ����Щ�������¸�ֵ
eval(['f=@(t,x,xd,u,ud,rd,p)',changeExp1(f_)]);
eval(['fu=@(t,x,xd,u,ud,rd,p)',changeExp1(fu_)]);
eval(['fud=@(t,x,xd,u,ud,rd,p)',changeExp1(fud_)]);
eval(['fx=@(t,x,xd,u,ud,rd,p)',changeExp1(fx_)]);
eval(['fxd=@(t,x,xd,u,ud,rd,p)',changeExp1(fxd_)]);
eval(['fr=@(t,x,xd,u,ud,rd,p)',changeExp1(fr_)]);

save(savemat_name);


function [x,xtheta,xtao,xtaou]=z_to_x(z,Dim);
    dimx=Dim.x;dimtheta=Dim.theta;dimtao=Dim.tao;
    x=z(1:dimx);s=dimx;
    shape=[dimx,dimtheta];
    xtheta=reshape(z(s+1:s+prod(shape)),shape);s=s+prod(shape);
    shape=[dimx,dimtao];
    xtao=reshape(z(s+1:s+prod(shape)),shape);s=s+prod(shape);
    shape=[dimx,1];
    xtaou=reshape(z(s+1:s+prod(shape)),shape);s=s+prod(shape);
function [xd,xthetad,xtaod,xtaoud,xdd]=zd_to_x(zd,Dim)
dimx=Dim.x;dimtheta=Dim.theta;dimtao=Dim.tao;
%with time delay tao
s=0;
xd=zd(s+1:s+dimx,1:dimtao);s=s+dimx;
xthetad=reshape(zd(s+1:s+dimx*dimtheta,1:dimtao),[dimx,dimtheta,dimtao]);s=s+dimx*dimtheta;
xtaod=reshape(zd(s+1:s+dimx*dimtao,1:dimtao),[dimx,dimtao,dimtao]);s=s+dimx*dimtao;
xtaoud=zd(s+1:s+dimx,1:dimtao);
%with time-delay tao+tao
xdd=zd(1:dimx,dimtao+1:dimtao+dimtao*dimtao);%%%
xdd=reshape(xdd',dimtao,dimtao,dimx);
xdd=permute(xdd,[3,1,2]);
function v=matrix1(symbol,n,m,k)
if nargin==2;
    v=sym(symbol,[n,1]);v=sym(v(:),'real');
elseif nargin==3;
    v=sym(symbol,[n,m]);v=sym(v,'real');
elseif nargin==4;
    for i=1:k;
        vi=matrix3d(symbol,n,m,k);
    end;
end
function x=vector(symbol,nx);
fsym=@(x,i)sym([x,num2str(i)],'real');
for i=1:nx;
    x(i)=fsym(symbol,i);
end;
x=x';
function x=matrix2d(symbol,n,m);
fsym2=@(x,i,j)sym(strcat(x,num2str(i),'_',num2str(j)),'real');
for i=1:n
    for j=1:m
        x(i,j)=fsym2(symbol,i,j);
    end;
end;
function x=matrix3d(symbol,n,m,k);
fsym3=@(x,i,j,l)sym(strcat(x,num2str(i),'_',num2str(j),'_',num2str(l)),'real');
for i=1:n;for j=1:m;for l=1:k;
            x(i,j,l)=fsym3(symbol,i,j,l);
        end;end;end;
function v=myjacobian(y,x)
ly=length(size(y));lx=numel(size(x));
if ly==2&size(y,2)==1
    if lx==2&size(x,2)==1;
        v=jacobian(y,x);
    elseif lx==2;
        v=jacobian2d(y,x);
    elseif lx==3;
        v=jacobian3d(y,x);
    end;
else
    'there is not appropriate method for this situation',
    wrong
end;
function v=jacobian3d(y,x);
m=size(x);
for s=1:length(y);
for i=1:m(1);
    for j=1:m(2);
        for k=1:m(3);
            v(s,i,j,k)=diff(y(s),x(i,j,k));
        end;
    end;
end;
end;
function v=jacobian_mod(y,x);
syms v;for i=1:length(x); if strcmp(class(x(i)),'sym');v(:,i)=diff(y,x(i));else v(:,i)=0;end;end;
function g=jacobian2d(f,xd);
for k=1:length(f);
for i=1:size(xd,1);
    for j=1:size(xd,2);
        g(k,i,j)=diff(f(k),xd(i,j));
    end;
end;
end;
function newstr=changeExp1(strMatrix)
 if isstruct(strMatrix)
     if numel(strMatrix)==0;
         newstr=strMatrix;
     else
        li=fieldnames(strMatrix);
        for i=1:numel(li);
            eval(['newstr.',li{i},'=','changeExp1(strMatrix.',li{i},');']);
        end;
     end;
elseif iscell(strMatrix)
    m=length(strMatrix);
    newstr={};
    for i=1:m;
        newstr{i}=changeExp1(strMatrix{i});
    end;
 else 
        s=strMatrix(:);
        if numel(s)==0;
            newstr='';
        else
            news=strcat('[',changeStr(char(s(1))));
            for i=2:numel(s);
                news=[news,'\n',changeStr(char(s(i)))];
            end;
            news=[news,']'];
            newstr=sprintf(news);%
        end;
end;
function newstr=changeStr(str)
str(str+0==32)='';%ȥ���ַ����еĿո�
n=length(str);
if n>0;
    start=[];finish=[];
    for i=2:n;
        if isnumber(str(i))&ischaracter(str(i-1)); 
           start=[start,i];
           if i<n;
             for j=i:n-1;
                 if j<n-1;
                    if isfinish(str(j+1));
                        finish=[finish,j];
                        break;
                    end;
                 else
                     finish=[finish,n];
                 end;
             end
           else
               finish=[finish,n];
          end;
       end;
    end;
    if length(finish)<length(start)
        str;
    end;

    move=0;newstr=str;
    for i=1:length(start);
        j=start(i);k=finish(i);
        varIndice=regexp(str(j:k),'_','split');%����'_'��ʶ���±�
        c=change(varIndice);
        a=start(i)+move;b=finish(i)+move;
        newstr=insert(newstr,a,b,c);
        move=move+(length(c)-b+a-1);
    end;
else
    newstr='';
end;
function newstr=insert(str,j,k,c);
newstr=[str(1:j-1),c,str(k+1:end)];
function c=change(varIndice)
dim=length(varIndice);
c=strcat('(',varIndice{1});
for i=2:dim;
    c=strcat(c,',',varIndice{i});
end;
c=strcat(c,')');
function flag=isnumber(c);%�Ƿ��Ǵ�0����ֵ
flag=(49<=c+0&c+0<=57);
function flag=isfinish(c);
flag=(48>c+0|c+0>57)&c+0~=95;
function flag=ischaracter(c);
flag=(c+0>=97&c+0<=122);

