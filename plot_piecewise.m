
function [h,tt,uu]=plot_piecewise(t,u,string)
%u是竖着的
%t的个数比u的个数要多1
% figure;
tt=linspace(min(t),max(t),1000);
uu=zeros(numel(tt),size(u,2));
n=numel(tt);
k=1;
for i=1:n;
    uu(i,:)=u(k,:);
    if (i+1<n & tt(i+1)>=t(k+1))
        k=k+1;
    end;
end;
h=[];
% if nargout==1;
    try
    h=plot(tt,uu,string,'linewidth',2);
    catch 
       h=plot(tt,uu,'linewidth',2);
    end;
% end;