function u1=interp_piecewise(tt,uu,t1)
k=dingwei(t1,tt);
if k==numel(tt);u1=uu(:,end);
else
    u1=uu(:,k);
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