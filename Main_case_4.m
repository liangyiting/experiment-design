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
model_id=1;
switch model_id;
    case 1
        load('model_1.mat');%model_1_symbol.mat
        xhistory=[2,0]';
        uhistory=[1;1];
        theta=[1,1,1,1];tao=1;taou=0;parameter=[theta,tao,taou];
        hx=[1 0;0 1];
        invQ=diag((1/0.1./[ 10 10]).^2);%�������
        uminvalue=[0.05;5];
        umaxvalue=[0.20;35];xmin=[1,0];xmax=[10,0.0000001];
        rdvalue=zeros(dimrd,1);
    case 2
        load('model_2.mat');
        xhistory=[3.3e-4,4e-3]';
        uhistory=[1e+5;4.4e+5];
        theta=[1,1];tao=1;taou=0;parameter=[theta,tao,taou];
        hx=[1 0;0 1];
        invQ=diag([1,1])/0.000003;%�������
        uminvalue=[0.6e+5;4.1e+5];
        umaxvalue=[1.5e+5;4.7e+5];
        rdvalue=[6e-4 9e-3]';
end
%%       
T=48;
dimy=size(hx,1);
Nsample=48;Nu=24;
tsample=linspace(0,T,Nsample+1);tsample=tsample(2:end);%������
measure_system=struct('tsample',tsample,'ysample',zeros(dimy,numel(tsample)),'hx',hx,'invQ',invQ);
tseq=linspace(0,T,Nu+1); 
% tt=linspace(0,T,200);dt=4*T/numel(tt);tt1=[];for i=1:numel(tseq);tt1=[tt1,linspace(max(0,tseq(i)-dt),min(T,tseq(i)+dt),10)];end;tt=unique(union(tt,tt1(tt1>=0&tt1<=T)));
tseq_rd=[0,tsample];rdseq_ref=rdvalue*ones(1,numel(tseq_rd)-1);rdfun_ref=@(t)interp_piecewise(tseq_rd,rdseq_ref,t);%�ο������ź�
rdseq=rdseq_ref;rdfun=@(t)interp_piecewise(tseq_rd,rdseq,t);
pestimate=parameter*1.0;pestimate(dimtheta+dimtao+1:end)=parameter(dimtheta+dimtao+1:end);%�������֪����ñ�ʶ��������ڹ����Ż�����Ϣ�Ż�
umin=uminvalue*ones(1,Nu);%umin=umin(:);
umax=umaxvalue*ones(1,Nu);%umax=umax(:);
%% ʵ���Ż����-4��ģʽ
if 0
    global criterion;
    criterion='A';
    urand=umin+rand(2,Nu).*(umax-umin);
    [~,~,~,li1,H1]=ObjForInfoOptim(urand,tseq,rdfun,parameter,xhistory);  
    fun=@(urand)ObjForInfoOptim(urand,tseq,rdfun,parameter,xhistory);
    ep=1e-6;[gdiff,df]=check_diff(fun,urand,1,ep);
end
global len_p;%�����ʶʵ������п��ǵĴ���ʶ�����ĳ���
global criterion;
criterion='A';mode=4;%mode=4����·��Լ������ʼ���Ż�ģʽģʽ %mode=3����·��Լ���Ż�


len_p=dimtheta;s=len_p;MaxIter=500;mode=4;%AN ������ʱ�Ͳ���
tic,[uoptinfo,fopt,u_kexing,xopt_a]=opt_Information_con(pestimate,rdfun,tseq,umin(:),umax(:),xhistory,uhistory,mode,MaxIter); toc
uopt_a=reshape(uoptinfo,2,Nu);
% tic,[uoptinfo1,fopt,u_kexing,xopt_a1]=opt_Information_con(pestimate,rdfun,tseq,umin(:),umax(:),xhistory,uhistory,mode,50); toc

len_p=dimtheta+dimtao;s=len_p;MaxIter=500;mode=4;% AD ����ʱ�Ͳ��� 
tic,[uoptinfo,fopt,u_kexing,xopt_b]=opt_Information_con(pestimate,rdfun,tseq,umin(:),umax(:),xhistory,uhistory,mode,MaxIter); toc
uopt_b=reshape(uoptinfo,2,Nu);


global x1max x2max epsilon1 epsilon2;
x1max=30;epsilon1=0.1;
x2max=10;epsilon2=0.1;
len_p=dimtheta+dimtao;s=len_p;
xd=[4.52,0]';MaxIter=500;mode=3; % AP ����ʱ�Ͳ�������������״̬Լ��
tic,[uoptinfo,fopt,u_kexing,xopt_b]=opt_Information_con(pestimate,rdfun,tseq,umin(:),umax(:),xd,uhistory,mode,MaxIter); toc
uopt_d=reshape(uoptinfo,2,Nu);
save('uopt_pathcon.mat','uopt_d');
tic,[~,~,~,li_d,H]=ObjForInfoOptim(uopt_d,tseq,rdfun,parameter,xd,uhistory);toc
figure;
for i=1:dimx;
    subplot(dimx,1,i);
    plot(li_d.t(li_d.t>=0),li_d.x(i,li_d.t>=0),'b-','linewidth',2);
    hold on;
    if i==1;plot([0,T],[x1max,x1max],'r-');ylabel('x_{1}[g/L]');axis([0,T,0,60])
    else plot([0,T],[x2max,x2max],'r-');ylabel('x_{2}[g/L]');xlabel('time[h]');axis([0,T,0,20])
    end
    hold off;
end
save('Main_case_4_0608.mat');
%% ���ؿ���ģ�� 4��ģʽ10��ģ���ʶ
%�����ʶ����
ubinary=umin;ubinary(1,2:2:end)=umax(1,2:2:end);ubinary(2,1:2:end)=umax(2,1:2:end);
pmin=max(0,parameter(1:dimtheta+dimtao)-8);pmax=parameter(1:dimtheta+dimtao)+8; %����ȫ���ǷǸ��ģ������ʹ���������
ptrue=parameter(1:dimtheta+dimtao);
udict={ubinary,uopt_a,uopt_b,uopt_d};
xdict={xhistory,xopt_a,xopt_b,xd};
cpdict={};
for kss=1:4;
    x=xdict{kss};
    u=udict{kss};
    [~,yout]=simulation(u,tseq,rdseq,tseq_rd,parameter,x,uhistory);
    cpopt=[];
    maxk=10;%10�����ؿ����ʶ����
    for k=1:maxk;
        noiseVar=sqrt(diag(inv(invQ)));%ֻ���ǲ���ص��ź�
        dimy=size(hx,1);ysample=yout;
        for i=1:dimy;ysample(i,:)=yout(i,:)+max(min(3*noiseVar(i),noiseVar(i)*randn(1,size(yout,2))),-3*noiseVar(i));end;
        MaxIter=100;
        popt=identify(ysample,u,rdseq_ref,s,pmin,pmax,parameter,tseq,tseq_rd,x,uhistory,MaxIter);
        cpopt=[cpopt;popt(:)'];
    end;
    cpdict{kss}=cpopt;
end
save('Main_case_4_0607_cpdict.mat');
%% table 3 ͳ�Ʋ���ƽ��ֵ�ͷ���
%ͳ������
[~,~,~,li,H]=ObjForInfoOptim(ubinary,tseq,rdfun,parameter,xhistory,uhistory);
[~,~,~,liopt_a,Hopt_a]=ObjForInfoOptim(uopt_a,tseq,rdfun,parameter,xopt_a,uhistory);
[~,~,~,liopt_b,Hopt_b]=ObjForInfoOptim(uopt_b,tseq,rdfun,parameter,xopt_b,uhistory);
[~,~,~,liopt_d,Hopt_d]=ObjForInfoOptim(uopt_d,tseq,rdfun,parameter,xd,uhistory);
fa=@(H)trace(inv(H(1:dimtheta+dimtao,1:dimtheta+dimtao))); %A  �о�
fe=@(H)sqrt(eig(inv(H(1:dimtheta+dimtao,1:dimtheta+dimtao))));% ��������
[fa(H),fa(Hopt_a),fa(Hopt_b),fa(Hopt_d)]
[fe(H),fe(Hopt_a),fe(Hopt_b),fe(Hopt_d)]

cstd=[];
mp=[];
for i=1:4;
    cpopt=cpdict{i};
    %stdi=sqrt(mean((cpopt-ones(maxk,1)*ptrue(:)').^2,1))%ͳ���������
    stdi=(std(cpopt,1));%ͳ�Ʊ�׼��
    cstd(i,:)=stdi;
    mp(i,:)=mean(cpopt,1);
end
cstd'
mp
figure;
for i=1:5;
    if i==5;subplot(3,2,1); 
    else subplot(3,2,i+2);
    end
    pi=[cpdict{1}(:,i),cpdict{2}(:,i),cpdict{3}(:,i)];
    plot(pi,'o-','markersize',5);hold on;plot([0,size(pi,1)],[1,1],'k--','linewidth',1.5);hold off;
    if i<5;title(['\vartheta_',num2str(i)])
    else title('\tau');
    end
    if i==5;legend('BS','AN','AD');end
end
%% figure 6
%��ͼ--����ģʽ�µõ�������״̬�켣ͼ figure 6
figure;
for i=1:dimx
    subplot(2,1,i);hold on;
    plot(li.t(li.t>=0),li.x(i,li.t>=0),'linewidth',2);
    plot(liopt_a.t(liopt_a.t>=0),liopt_a.x(i,liopt_a.t>=0),'r-','linewidth',2)
    plot(liopt_b.t(liopt_b.t>=0),liopt_b.x(i,liopt_b.t>=0),'k-','linewidth',2)
    plot(liopt_d.t(liopt_d.t>=0),liopt_d.x(i,liopt_d.t>=0),'g-','linewidth',2)
    if i==1;ylabel('x_{1}(g/L)');
    else;ylabel('x_{2}(g/L)');
    end
    hold on;
    if i==1;plot([0,T],[x1max,x1max],'k--');ylabel('x_{1}[g/L]');
    else plot([0,T],[x2max,x2max],'k--');ylabel('x_{2}[g/L]');xlabel('time[h]')
    end
    hold off;
    xlabel('Time(h)');legend('BS','AN','AD','AP'); 
end

%% figure 5
figure;%4��ģʽ�Ŀ��ƹ켣ͼ
subplot(1,2,1);plot_piecewise(tseq,ubinary(1,:)','b-');ylabel('u[h^{-1}]');
axis([0,T,umin(1),umax(1)*1.1])
subplot(1,2,2);plot_piecewise(tseq,ubinary(2,:)','b-');ylabel('v[g/L]');
axis([0,T,umin(2),umax(2)*1.1])
% title('BS');
figure;
subplot(1,2,1);plot_piecewise(tseq,uopt_a(1,:)','b-');ylabel('u[h^{-1}]');
hold on;plot([0,T],[umaxvalue(1),umaxvalue(1)],'r--');hold off;
axis([0,T,umin(1),umax(1)*1.1])
subplot(1,2,2);plot_piecewise(tseq,uopt_a(2,:)','b-');ylabel('v[g/L]');
hold on;plot([0,T],[umaxvalue(2),umaxvalue(2)],'r--');hold off;
axis([0,T,umin(2),umax(2)*1.1])
% title('AN')
figure;
subplot(1,2,1);plot_piecewise(tseq,uopt_b(1,:)','b-');ylabel('u[h^{-1}]');
hold on;plot([0,T],[umaxvalue(1),umaxvalue(1)],'r--');hold off;
axis([0,T,umin(1),umax(1)*1.1])
subplot(1,2,2);plot_piecewise(tseq,uopt_b(2,:)','b-');ylabel('v[g/L]');
hold on;plot([0,T],[umaxvalue(2),umaxvalue(2)],'r--');hold off;
axis([0,T,umin(2),umax(2)*1.1])
% title('AD')
figure;
subplot(1,2,1);plot_piecewise(tseq,uopt_d(1,:)','b-');ylabel('u[h^{-1}]');
hold on;plot([0,T],[umaxvalue(1),umaxvalue(1)],'r--');hold off;
axis([0,T,umin(1),umax(1)*1.1])
subplot(1,2,2);plot_piecewise(tseq,uopt_d(2,:)','b-');ylabel('v[g/L]');
hold on;plot([0,T],[umaxvalue(2),umaxvalue(2)],'r--');hold off;
axis([0,T,umin(2),umax(2)*1.1])

