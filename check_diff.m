function [gdiff,df]=check_diff(fun,useq,s,ep)
[f0,df]=fun(useq);
dimu=size(useq,1);
gdiff=useq*0;
try
    ep;
catch
    ep=1e-6;
end
try 
    s;
catch
    s=3;
end;
for i=1:dimu*s%numel(useq);
    ui=useq;ui(i)=useq(i)+ep;gdiff(i)=(fun(ui)-f0)/ep;
end;
format long
try;
[gdiff(:,1:s),df(:,1:s)]
end;

