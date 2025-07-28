function [ishk] = SENSOR(u)

% Set up sensor
il=size(u,1);

psi=zeros(il,1);

% WENO Smoothness sensor on pressure
a=1/4;
b=13/12;
nshrp=2;
nbuf=3;
eps=0.000005; 
eps=1e-5;    % <-- WENO type sensor threshold
% eps=-1000; % <-- use WENO everywhere

ii=1;
phi = (a*1/4*(u(ii+1) - u(il-1))^2 +  b*(u(ii+1) - 2*u(ii)+ u(il-1))^2 )^nshrp; %I,I-1,I+1 stencil

% num = u(ii+1) - 2*u(ii) + u(il-1);
% den = u(ii+1) + 2*u(ii) + u(il-1);
% phi = abs(num)/abs(den);

if phi>eps
  % psi(ii) = 1; % Include buffer points
  ilow=max(ii-nbuf,1);
  imax=min(ii+nbuf,il);
  psi(ilow:imax) = 1; % Include buffer points
end



for ii=2:il-1
  
 phi = (a*1/4*(u(ii+1) - u(ii-1))^2 +  b*(u(ii+1) - 2*u(ii)+ u(ii-1))^2 )^nshrp; %I,I-1,I+1 stencil

 % num = u(ii+1) - 2*u(ii) + u(ii-1);
 % den = u(ii+1) + 2*u(ii) + u(ii-1);
 % phi = abs(num)/abs(den);

 if phi>eps
   ilow=max(ii-nbuf,1);
   imax=min(ii+nbuf,il);
   psi(ilow:imax) = 1; % Include buffer points
 end

end

ii=il;
phi = (a*1/4*(u(2) - u(il-1))^2 +  b*(u(2) - 2*u(ii)+ u(il-1))^2 )^nshrp; %I,I-1,I+1 stencil

% num = u(2) - 2*u(ii) + u(ii-1);
% den = u(2) + 2*u(ii) + u(ii-1);
% phi = abs(num)/abs(den);

if phi>eps
  ilow=max(ii-nbuf,1);
  imax=min(ii+nbuf,il);
  psi(ilow:imax) = 1; % Include buffer points
  % psi(ii) = 1; % Include buffer points
end

for ii=1:il-1
  psihat(ii)=psi(ii+1)-psi(ii);
end

ishk=find(psi==1);


end