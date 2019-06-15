function [ysample,yout]=simulation(useq,tseq,rdseq,tseq_rd,parameter,xhistory,uhistory)
global measure_system;
[~,~,tt,xx]=ObjForIdentify_1(parameter,rdseq,tseq_rd,useq,tseq,xhistory,uhistory, measure_system);
hx=measure_system.hx;tsample=measure_system.tsample;invQ=measure_system.invQ;
yy=hx*xx(1:2,:);yout=interp1(tt,yy',tsample)';if size(yout,2)==1;yout=yout';end;
noiseVar=sqrt(diag(inv(invQ)));%只考虑不相关的信号
dimy=size(hx,1);ysample=yout;
function [J,dJp,tt,xx]=ObjForIdentify_1(parameter,rdseq,tseq_rd,useq,tseq,xhistory,uhistory, measure_system)
if nargout>1;
    [J,dJp,tt,xx]=ObjForIdentify(rdseq,tseq_rd,parameter,useq,tseq,xhistory,uhistory,measure_system);
else
    J=ObjForIdentify(rdseq,tseq_rd,parameter,useq,tseq,xhistory,uhistory, measure_system);
end 
