function [IC] = FDICMKR(xv,a0,a1)
% -- returns FD initial condition and 1D grid
%%% Point Value Grid
% dx=(xf-x0)/(nx-1);
% xv=x0:dx:xf;
% xc=x0+dx/2:dx:xf-dx/2;
% xc=xc';

nx=length(xv);

%%% Creating IC
a0 = a0; %0.5;
a1 = a1; %1;
IC = zeros(nx,1);
IC = a0 + a1*sin(xv);


end