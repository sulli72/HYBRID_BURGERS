function [RHS] = XWENOFLUX(u,xv)
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

%%% Split  reconstructed flux vectors
RGHTFLX=zeros(il,1);
LEFTFLX=zeros(il,1);

%%% Obtain flux vector
fv=U2FLX(u);

%%% Split the Flux (Locally or Globally?)
[fpls, fmns] = LAXFRD_SPLIT(fv,u);


%%% Visualize split flux if needed
if iplt==1
  figure(7)
  plot(xv,fv,xv,fpls,'--',xv,fmns,'--',xv,fpls+fmns,'o','LineWidth',2)
  grid on
  legend('$f(u)$','$f^+(u)$','$f^-(u)$','$f^+ + f^-$','Location','best')
  xlabel('$x$')
  ylabel('$f(u)$')
  db=1;
end

% WENO Procedure for f+ fluxes
% uses i-2,i-1,i,i+1,i+2 points

% Linear weights for O(dx^3) stencil
d0=3/10; %3/10;
d1=3/5;
d2=1/10; %1/10;

eps=1E-8; % div. by 0 fix

lbb=1;
lim1=[il-1 il];
lim2=[il-2 il-1];

rbb=1;
rip1=[1 2];
rip2=[2 3];

fm2vec=circshift(fpls,2);
fm1vec=circshift(fpls,1);
fp1vec=circshift(fpls,-1);
fp2vec=circshift(fpls,-2);

for ii=1:il

  if ii<=2
    fm2=fpls(lim2(lbb));
    fm1=fpls(lim1(lbb));
    fc0=fpls(ii);
    fp1=fpls(ii+1);
    fp2=fpls(ii+2);

    lbb=lbb+1;
  elseif ii>=il-1
    fm2=fpls(ii-2);
    fm1=fpls(ii-1);
    fc0=fpls(ii);
    fp1=fpls(rip1(rbb));
    fp2=fpls(rip2(rbb));

    rbb=rbb+1;
  else
    fm2=fpls(ii-2);
    fm1=fpls(ii-1);
    fc0=fpls(ii);
    fp1=fpls(ii+1);
    fp2=fpls(ii+2);
  end

  % Circular shift approach
  %     fm2=fm2vec(ii);
  %     fm1=fm1vec(ii);
  %     fc0=fpls(ii);
  %     fp1=fp1vec(ii);
  %     fp2=fp2vec(ii);

  % Weights for ENO Stencils
  c0j = [1/3 5/6 -1/6]; % I,I+1,I+2 stencil (r=0 stencil)
  c1j = [-1/6 5/6 1/3]; % I-1,I,I+1 stencil (r=1 stencil)
  c2j = [1/3 -7/6 11/6]; % I-2,I-1,I stencil (r=2 stencil)

  S0 = [fc0 fp1 fp2];
  S1 = [fm1 fc0 fp1];
  S2 = [fm2 fm1 fc0];

  fr0 = c0j*S0';
  fr1 = c1j*S1';
  fr2 = c2j*S2';

  % Compute smoothness indicator polynomials
  B0=13/12*(fc0 - 2*fp1 + fp2)^2 + 1/4*(3*fc0 - 4*fp1 + fp2)^2; %I,I-1,I-2 stencil (measures 1st, 2nd derivs using FD of stencil cell avgs)
  B1=13/12*(fm1 - 2*fc0 + fp1)^2 + 1/4*(fm1 - fp1)^2; %I,I-1,I+1 stencil
  B2=13/12*(fm2 - 2*fm1 + fc0)^2 + 1/4*(fm2 - 4*fm1 + 3*fc0)^2; %I,I+1,I+2 stencil

  % Compute alpha values (needed for non-linear weights)
  a0 = d0/(eps+B0)^2;
  a1 = d1/(eps+B1)^2;
  a2 = d2/(eps+B2)^2;
  as = a0 + a1 + a2;

  % Compute Non-Linear Weights
  w0 = a0/as;
  w1 = a1/as;
  w2 = a2/as;

  ws=w0+w1+w2;
  if abs(1-ws) > 0.00001
  fprintf(1,'Bad Weight Combo, Wsum = %f\n',ws)
  end
  RGHTFLX(ii) = w0*fr0 + w1*fr1 + w2*fr2;

end


if iplt==1
  %%% Visualize RECONSTRUCTED!!! split flux if needed
  figure(8)
  plot(xv,fv,xv,fpls,'--',(xv+0.5*dx),RGHTFLX,'o','LineWidth',2)
  grid on
  legend('$f(u)$','$f^+(u)$','$f^+_{WENO}$','Location','best')
  xlabel('$x$')
  ylabel('$f(u)$')
  db=1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% ------- FLUX FOR LEFT RUNNING WAVES ----- %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% WENO using appropriate stencils for the flux from the left running wave, f-

% Linear weights for O(dx^3) stencil
d0=3/10; % ENSURE LINEAR WEIGHTS ARE CORRECT FOR f- TERMS
d1=3/5;
d2=1/10;

% Weights for ENO Stencils
cn1j = [11/6 -7/6 1/3]; % I+1,I+2,I+3 stencil
c0j = [1/3 5/6 -1/6]; % I-2,I-1,I stencil
c1j = [-1/6 5/6 1/3]; % I-1,I,I+1 stencil %%% THIS IS WRONG IN UMNS CODE??
c2j = [1/3 -7/6 11/6]; % proper weights for I+1,I+2,I+3 stencil?

% Boundary point stencils
ip1=[il-1 1 2];
ip2=[1 2 3];
ip3=[2 3 4];
bb=1;

fm1vec=circshift(fmns,1);
fp1vec=circshift(fmns,-1);
fp2vec=circshift(fmns,-2);
fp3vec=circshift(fmns,-3);

for ii=1:il

    if ii==1
  
      fm1=fmns(il-1);
      fc0=fmns(ii);
      fp1=fmns(ii+1);
      fp2=fmns(ii+2);
      fp3=fmns(ii+3);
  
      lbb=lbb+1;
    elseif ii>=il-2
  
      fm1=fmns(ii-1);
      fc0=fmns(ii);
      fp1=fmns(ip1(bb));
      fp2=fmns(ip2(bb));
      fp3=fmns(ip3(bb));
      bb=bb+1;
    else
      fm1=fmns(ii-1);
      fc0=fmns(ii);
      fp1=fmns(ii+1);
      fp2=fmns(ii+2);
      fp3=fmns(ii+3);
    end

  %   % Circular shift approach
%   fm1=fm1vec(ii);
%   fc0=fmns(ii);
%   fp1=fp1vec(ii);
%   fp2=fp2vec(ii);
%   fp3=fp3vec(ii);

  % Weights for ENO Stencils
  c0j = [1/3 5/6 -1/6]; % I+1,I,I-1 stencil (r=0 stencil)
  c1j = [-1/6 5/6 1/3]; % I+2,I+1,I stencil (r=1 stencil)
  c2j = [1/3 -7/6 11/6]; % I+3,I+2,I+1 stencil (r=2 stencil)

  S0 = [fp1 fc0 fm1];
  S1 = [fp2 fp1 fc0];
  S2 = [fp3 fp2 fp1];

%   fvup = [fm1 fc0 fp1];
%   fvcn = [fc0 fp1 fp2];
%   fvdn = [fp1 fp2 fp3];

  fr0 = c0j*S0';
  fr1 = c1j*S1';
  fr2 = c2j*S2';
  
%   fintm = cr0*fvup';
%   fintc = cr1*fvcn';
%   fintp = cr2*fvdn';

  % Compute smoothness indicator polynomials
  B0=13/12*(fp1 - 2*fc0 + fm1)^2 + 1/4*(3*fp1 - 4*fc0 + fm1)^2; %I-1,I,I+1 stencil (measures 1st, 2nd derivs using FD of stencil cell avgs)
  B1=13/12*(fp2 - 2*fp1 + fc0)^2 + 1/4*(fp2 - fc0)^2; %I,I+1,I+2 stencil
  B2=13/12*(fp3 - 2*fp2 + fp1)^2 + 1/4*(fp3 -4*fp2 + 3*fp1)^2; %I+1,I+2,I+3 stencil


  % Compute alpha values (needed for non-linear weights)
  a0 = d0/(eps+B0)^2;
  a1 = d1/(eps+B1)^2;
  a2 = d2/(eps+B2)^2;
  as = a0 + a1 + a2;

  % Compute Non-Linear Weights
  w0 = a0/as;
  w1 = a1/as;
  w2 = a2/as;

  ws=w0+w1+w2;
  if abs(1-ws) > 0.00001
  fprintf(1,'Bad Weight Combo, Wsum = %f\n',ws)
  end

  LEFTFLX(ii) = w0*fr0 + w1*fr1 + w2*fr2;

end
% LEFTFLX=circshift(LEFTFLX,1);
FLUXFACE=RGHTFLX+LEFTFLX;

%%% Visualize RECONSTRUCTED!!! split flux if needed
if iplt==1
  figure(9)
  plot(xv,fv,xv,fmns,'--',(xv+0.5*dx),LEFTFLX,'o','LineWidth',2)
  grid on
  legend('$f(u)$','$f^+(u)$','$f^-_{WENO}$','Location','best')
  xlabel('$x$')
  ylabel('$f(u)$')
  db=1;


  figure(10)
  plot(xv,fv,xv,fmns,'--',xv,fpls,'--',(xv+0.5*dx),FLUXFACE,'o','LineWidth',2)
  grid on
  legend('$f(u)$','$f^-(u)$','$f^+$','$\hat{f}_{j+{1/2}}$','Location','best')
  xlabel('$x$')
  ylabel('$f(u)$')
  db=1;
end

FP1=FLUXFACE;
FM1=circshift(FP1,1);
FM1(1)=FM1(il); % PERIODIC BC HARD FIX
RHS=FP1-FM1;
db=1;
end