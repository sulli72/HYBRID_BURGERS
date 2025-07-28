function [linf,l1,l2,uex] = FDERRCOMP(xv,tend,unum,a0,a1)
% Function that returns error between exact and numerical solutions for
% NONLINEAR PROBLEMS.
% for linear, scalar equations might be easier to just use method of
% characteristics...

% a0=0.5;
% a1=1;

nx=length(xv);
dx=xv(2)-xv(1);
uex=zeros(nx,1);
t=tend;

np=length(xv);
uex=zeros(np,1);


for i=1:nx
  x=xv(i);


  tol=1E-6;
  maxiters=100;
  iters=0;
  phinew=x-x/100;
  res=1E8;

  while res>tol && iters<maxiters
    phi=phinew;
    g = x-phi-(a0 + a1*sin(phi))*t;
    gpr = -1 - sign(a1)*cos(phi)*t;

    phinew = phi - sign(a1)*g/gpr;

    res=abs(phinew-phi);
    iters=iters+1;

  end

  if iters>=maxiters
    fprintf(1,'Exceeded Newton iteration limit\n')
  end


  root=phinew;
  U = a0 + a1*sin(root);
  uex(i) = U;

end

%% Simpsons Method for Cell Average Integration Step
% jj=1;
% for ii=2:2:nx
%   IU = dx/3*(uex(ii-1) + 4*uex(ii) + uex(ii+1));
%   uex(jj)=IU/(2*dx);
%   jj=jj+1;
% end




% figure(6)
% plot(xv,uex,'-')
% grid on
% db=1;

% Norm Calculation
diff=uex-unum;
linf=max(abs(diff));
l1=dx*sum(abs(diff));
l2=sqrt(dx)*sqrt(sum(diff.^2));

end