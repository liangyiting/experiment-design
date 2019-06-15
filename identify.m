function popt=identify(ysample,useq,rdseq_ref,s,pmin,pmax,parameter,tseq,tseq_rd,xhistory,uhistory,MaxIter)
global measure_system;
measure_system.ysample=ysample;
 p0=pmin(:)+rand(numel(pmax),1).*(pmax(:)-pmin(:));p0=p0(1:s);p0=p0(:);
 p0=parameter(1:s);p0=p0(:);
options=optimset('display','iter','algorithm','interior-point','GradObj','on');
try options.MaxIter=MaxIter;
catch options.MaxIter=100;
end
fun=@(p)ObjForIdentify(rdseq_ref,tseq_rd,[p(:)',parameter(s+1:end)],useq,tseq,xhistory,uhistory,measure_system);
if 0; p1=4*p0.*rand(size(p0));[a,b]=check_diff(fun,p1,1,1e-6); end
tic,[popt,fopt]=fmincon(fun,p0(:),[],[],[],[],pmin(:),pmax(:),[],options);toc
function [gdiff,df]=check_diff(fun,useq,s,ep)
[f0,df]=fun(useq);
dimu=size(useq,1);
gdiff=useq*0;
try;ep;catch;ep=1e-6;end
try; s;catch;s=3;end;
for i=1:dimu*s%numel(useq);
    ui=useq;ui(i)=useq(i)+ep;gdiff(i)=(fun(ui)-f0)/ep;
end;
format long
try;
[gdiff(:,1:s),df(:,1:s)]
end;
