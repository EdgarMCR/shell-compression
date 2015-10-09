function [r]=errorNorm(x, xnew, L)
%Returns the error given an exact solution and an approximate solution as
%well as the length L of the domain
    n=length(x);
    temp=0.0;
    for ii=1:n
        temp=temp+(x(ii)-xnew(ii))^2*(L/(n-1));
    end
    temp=temp/L;
    r=sqrt(temp);
end