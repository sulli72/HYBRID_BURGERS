function [RHS] = CENTRAL6(u,xv)
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
RHS(ii) = (-1*fv(il-3) + 9*fv(il-2) - 45*fv(il-1) + 45*fv(ii+1) - 9*fv(ii+2) + 1*fv(ii+3)) /(60);

ii=2;
RHS(ii) = (-1*fv(il-2) + 9*fv(il-1) - 45*fv(ii-1) + 45*fv(ii+1) - 9*fv(ii+2) + 1*fv(ii+3)) /(60);

ii=3;
RHS(ii) = (-1*fv(il-1) + 9*fv(ii-2) - 45*fv(ii-1) + 45*fv(ii+1) - 9*fv(ii+2) + 1*fv(ii+3)) /(60);

for ii=4:il-3
  RHS(ii) = (-1*fv(ii-3) + 9*fv(ii-2) - 45*fv(ii-1) + 45*fv(ii+1) - 9*fv(ii+2) + 1*fv(ii+3)) /(60);
end

ii=il-2;
RHS(ii) = (-1*fv(ii-3) + 9*fv(ii-2) - 45*fv(ii-1) + 45*fv(ii+1) - 9*fv(ii+2) + 1*fv(2)) /(60);

ii=il-1;
RHS(ii) = (-1*fv(ii-3) + 9*fv(ii-2) - 45*fv(ii-1) + 45*fv(ii+1) - 9*fv(2) + 1*fv(3)) /(60);

ii=il;
RHS(ii) = (-1*fv(ii-3) + 9*fv(ii-2) - 45*fv(ii-1) + 45*fv(2) - 9*fv(3) + 1*fv(4)) /(60);


end

