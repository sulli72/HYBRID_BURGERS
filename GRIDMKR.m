function [xv] = GRIDMKR(xf,x0,il,grdopt)
% -- returns 1D grid w/ desired spacing

% Copmute baseline dx
dx=(xf-x0)/(il-1);

if grdopt == 1
  % fprintf(1,'Doing constant spacing mesh\n')

  % Constant spacing grid for convergence checks
  dx=(xf-x0)/(il-1);
  xv=x0:dx:xf;

elseif grdopt == 2
  % fprintf(1,'Doing fine near shock spacing mesh\n')


  % Get sinusoidal mesh spacing since periodic problem
  Lx=xf-x0;
  DXAMP = 0.75; % <-- aggressiveness of mesh spacing change
  xshk=-pi/2;
  ishft=2*(xshk/Lx)*2*pi; % <-- shift smallest spacing to near shock
  
  xmy=x0;
  xv(1) = x0;
  for i=1:il-1
    xmy=xv(i);
    dxmy = dx*(1 + DXAMP*sin(2*pi*i/il - ishft));
    xv(i+1) = xmy+dxmy;
    dxv(i)=dxmy;
  end

  % Rescale so refinement ends exactly at xf
  xv = x0 + ( (xv-x0)./(xv(end) - x0) ).*(xf-x0);

elseif grdopt == 3
  % fprintf(1,'Doing coarse near shock spacing mesh\n')

    % Get sinusoidal mesh spacing since periodic problem
  Lx=xf-x0;
  DXAMP = 0.75; % <-- aggressiveness of mesh spacing change
  xshk=-pi/2;
  ishft=0; % <-- keep smallest spacing to away from shock
  
  xmy=x0;
  xv(1) = x0;
  for i=1:il-1
    xmy=xv(i);
    dxmy = dx*(1 + DXAMP*sin(2*pi*i/il - ishft));
    xv(i+1) = xmy+dxmy;
    dxv(i)=dxmy;
  end

  % Rescale so refinement ends exactly at xf
  xv = x0 + ( (xv-x0)./(xv(end) - x0) ).*(xf-x0);

elseif grdopt == 4
  % fprintf(1,'Doing ''random'' spacing mesh\n')

  % PERTURB MESH BY CERTAIN FACTOR 
  % factor chosen to maintain some degree of mesh continuity
  amp=0.4;

  % Constant spacing grid for convergence checks
  dx=(xf-x0)/(il-1);
  xv=x0:dx:xf;

  for i=2:il-1
    xrn = randn(1);
    sgn = sign(xrn);
    xv(i) = xv(i) + sgn*amp*dx;

  end

else
  fprintf(1,'INVALID MESHING OPTION\n')
end


end