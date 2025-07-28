function [fpls, fmns] = LAXFRD_SPLIT(fv,u)

nx=length(fv);
fpls=zeros(nx,1);
fmns=zeros(nx,1);

% Simple flux splitting for Burgers flux function
% f+ = max(0,u)*u/2
% f- = min(0,u)*u/2
alpha=max(abs(u));
for ii=1:nx

%  alpha=max(u);
 % Local flux splitting approach (at one grid point)
%  fpls(ii) = 0.5*u(ii)*max(0,u(ii));
%  fmns(ii) = 0.5*u(ii)*min(0,u(ii));

 % Global flux splitting approach (over whole domain)
 fpls(ii) = 0.5*(fv(ii) + alpha*u(ii));
 fmns(ii) = 0.5*(fv(ii) - alpha*u(ii));

end


end