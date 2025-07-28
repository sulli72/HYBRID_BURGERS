function [xix] = METRIC(xv)
% Compute 1D metrics for solving equation on stretched mesh

% Compute Metrics using 6th order central scheme

% use dxi = const = 1;

% PERIODICITY FIX - ADD/SUBTRACT LX AS NECESSARY
x0=xv(1);
xf=xv(end);
il=length(xv);

ii=1;
Lx = xf-x0;

xxi(ii) = (-1*( xv(il-3) - Lx ) + 9*( xv(il-2) - Lx ) - 45*( xv(il-1) - Lx ) + 45*xv(ii+1) - 9*xv(ii+2) + 1*xv(ii+3)) /(60);

ii=2;
xxi(ii) = (-1*( xv(il-2) - Lx ) + 9*( xv(il-1) - Lx ) - 45*( xv(ii-1)  ) + 45*xv(ii+1) - 9*xv(ii+2) + 1*xv(ii+3)) /(60);

ii=3;
xxi(ii) = (-1*( xv(il-1) - Lx ) + 9*( xv(ii-2)  ) - 45*( xv(ii-1)  ) + 45*xv(ii+1) - 9*xv(ii+2) + 1*xv(ii+3)) /(60);

for ii=4:il-3
  xxi(ii) = (-1*xv(ii-3) + 9*xv(ii-2) - 45*xv(ii-1) + 45*xv(ii+1) - 9*xv(ii+2) + 1*xv(ii+3)) /(60);
end

ii=il-2;
xxi(ii) = (-1*xv(ii-3) + 9*xv(ii-2) - 45*xv(ii-1) + 45*xv(ii+1) - 9*xv(ii+2) + 1*( Lx + xv(2) )) /(60);

ii=il-1;
xxi(ii) = (-1*xv(ii-3) + 9*xv(ii-2) - 45*xv(ii-1) + 45*xv(ii+1) - 9*( Lx + xv(2) ) + 1*( Lx + xv(3)) ) /(60);

ii=il;
xxi(ii) = (-1*xv(ii-3) + 9*xv(ii-2) - 45*xv(ii-1) + 45*( Lx + xv(2) ) - 9*( Lx + xv(3)) + 1*( Lx + xv(4) )) /(60);


xix = 1./xxi;

end