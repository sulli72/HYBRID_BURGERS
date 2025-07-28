function [RHS] = HYBRIDMASK(utmp,ishk,xv,n)

[RHS] = CENTRAL6(utmp,xv);

% inds=find(psi==1);
il=length(utmp);

% ADD GHOST POINTS FOR PERIODIC BOUNDARY CONDTIONS
ughost=[utmp(il-3:il-1); utmp; utmp(2:4)  ];

for ii=1:length(ishk)
  imy=ishk(ii);

  % Include enough points for WENO stencil to operate cleanly
  ilow=imy-3;
  iupp=imy+3;
  [RHS(imy)] = XWENOMASK(ughost(ilow+3:iupp+3));
end

if n==20 || n==100
  teleflux=sum(RHS(1:il-1));
  fprintf(1,'Telescope flux sum = %f\n',teleflux)
  db=1;
end



end