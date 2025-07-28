function [RHS] = UPWIND1(u,xv)
% Inputs : flux vector, grid vector
% Outputs: flux at the interfaces from the right:
% i.e. :

% j-1/2     j      j+1/2
%   |                |
%   |-------0--------| <--- FLUX CALCULATED HERE AT X^+_j+1/2
%   |                |
% This function computes the flux at the j+1/2 interface from the RIGHT
% (F^+_j+1/2)

%%% Plotting for vis/debug
iplt=0;

%%% Global constants
il=length(xv);
dx=xv(2)-xv(1);

% Allocate variable
RHS=zeros(il,1);

%%% Obtain flux vector
fv=U2FLX(u);

%%% Compute DFDX

ii=1;
RHS(ii) = (fv(ii) - fv(il-1) )/1;

for ii=2:il
  RHS(ii) = (fv(ii) - fv(ii-1) )/1;
end


end

