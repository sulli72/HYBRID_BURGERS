function PHIDERIV = COMPACTDIFF6(u,xv)
% 6th Order Compact FD scheme
il=length(u);

PHI=U2FLX(u);

dxi=1;

% Scheme coeffs
% C6 interior
alpha=1/3; a=14/9; b=1/9; 

% C4 i=1, i=il
alpha1=3;
a1=-17/6; b1=3/2; c1=3/2;
d1=-1/6; e1=0;

% AC5 i=2,i=il-1
alpha21=3/14; alpha22=3/14;
a2=-19/28; b2=-5/42; c2=6/7;
d2=-1/14; e2=1/84;

% Tri-diag matrix vectors
ap1=zeros(il-1,1); % super diag
ac0=ones(il,1);   % main diag -- hard fixed for dirichlet BCs
am1=zeros(il-1,1); % sub diag
kv=zeros(il,1);    % RHS vector

% i=1 point - C4 scheme
ii=1;
ap1(ii)=alpha1;
ac0(ii)=1;
kv(ii) = 1/dxi * (a1*PHI(ii) + b1*PHI(ii+1) + c1*PHI(ii+2) + d1*PHI(ii+3) + e1*PHI(ii+4));

% i=2 point -  AC5 scheme
ii=2;
ap1(ii)=alpha22;
ac0(ii)=1;
am1(ii-1)=alpha21;
kv(ii) = 1/dxi * (a2*PHI(ii-1) + b2*PHI(ii) + c2*PHI(ii+1) + d2*PHI(ii+2) + e2*PHI(ii+3));

% Build tri-diag vectors
for ii=3:il-2
  ap1(ii)=alpha;
  ac0(ii)=1;
  am1(ii-1)=alpha;
  kv(ii)=a*(PHI(ii+1)-PHI(ii-1))/(2*dxi) + b*(PHI(ii+2)-PHI(ii-2))/(4*dxi);
end


% i=il-1 point -  AC5 scheme
ii=il-1;
ap1(ii)=alpha21;
ac0(ii)=1;
am1(ii-1)=alpha22;
kv(ii) = -1/dxi * (a2*PHI(ii+1) + b2*PHI(ii) + c2*PHI(ii-1) + d2*PHI(ii-2) + e2*PHI(ii-3));


% i=il point - C4 scheme
ii=il;
am1(ii-1)=alpha1;
ac0(ii)=1;
kv(ii) = -1/dxi * (a1*PHI(ii) + b1*PHI(ii-1) + c1*PHI(ii-2) + d1*PHI(ii-3) + e1*PHI(ii-4));


PHIDERIV = THOMAS(ac0,ap1,am1,kv);
db=1;

end