function [r]=residueNorm(r, L)
%Returns the error given an exact solution and an approximate solution as
%well as the length L of the domain
    n=length(r);
    temp=0.0;
    for ii=1:n
        temp=temp+(r(ii))^2*(L/(n-1));
    end
    temp=temp/L;
    r=sqrt(temp);
end