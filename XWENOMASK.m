function [RHSLOC] = XWENOMASK(u)
% WENO Operation on small set of input variables
% Inputs : flux vector, grid vector
% Outputs: flux at the interfaces from the right:
% i.e. :

% j-1/2     j      j+1/2
%   |                |
%   |-------0--------| <--- FLUX CALCULATED HERE AT X^+_j+1/2
%   |                |
% This function computes the flux at the j+1/2 interface from the RIGHT
% (F^+_j+1/2)

%%% Obtain flux vector
fv=U2FLX(u);

%%% Split the Flux (Locally or Globally?)
[fpls, fmns] = LAXFRD_SPLIT(fv,u);


%% WENO F i+1/2 +/- STENCIL ON i-2 i-1 i i+1 i+2 i+3

% Linear weights for O(dx^3) stencil
d0=3/10; %3/10;
d1=3/5;
d2=1/10; %1/10;

eps=1E-8; % div. by 0 fix

fm2=fpls(2);
fm1=fpls(3);
fc0=fpls(4);
fp1=fpls(5);
fp2=fpls(6);

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
RGHTFLX = w0*fr0 + w1*fr1 + w2*fr2;




% if iplt==1
%   %%% Visualize RECONSTRUCTED!!! split flux if needed
%   figure(8)
%   plot(xv,fv,xv,fpls,'--',(xv+0.5*dx),RGHTFLX,'o','LineWidth',2)
%   grid on
%   legend('$f(u)$','$f^+(u)$','$f^+_{WENO}$','Location','best')
%   xlabel('$x$')
%   ylabel('$f(u)$')
%   db=1;
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% ------- FLUX FOR LEFT RUNNING WAVES ----- %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% WENO using appropriate stencils for the flux from the left running wave, f-

% Linear weights for O(dx^3) stencil
d0=3/10; % ENSURE LINEAR WEIGHTS ARE CORRECT FOR f- TERMS
d1=3/5;
d2=1/10;

fm1=fmns(3);
fc0=fmns(4);
fp1=fmns(5);
fp2=fmns(6);
fp3=fmns(7);


% Weights for ENO Stencils
c0j = [1/3 5/6 -1/6]; % I+1,I,I-1 stencil (r=0 stencil)
c1j = [-1/6 5/6 1/3]; % I+2,I+1,I stencil (r=1 stencil)
c2j = [1/3 -7/6 11/6]; % I+3,I+2,I+1 stencil (r=2 stencil)

S0 = [fp1 fc0 fm1];
S1 = [fp2 fp1 fc0];
S2 = [fp3 fp2 fp1];



fr0 = c0j*S0';
fr1 = c1j*S1';
fr2 = c2j*S2';

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

LEFTFLX = w0*fr0 + w1*fr1 + w2*fr2;


FP1=RGHTFLX+LEFTFLX;


%% WENO F i-1/2 +/- STENCIL ON i-3 i-2 i-1 i i+1 i+2

fm2=fpls(1);
fm1=fpls(2);
fc0=fpls(3);
fp1=fpls(4);
fp2=fpls(5);

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
RGHTFLX = w0*fr0 + w1*fr1 + w2*fr2;




% if iplt==1
%   %%% Visualize RECONSTRUCTED!!! split flux if needed
%   figure(8)
%   plot(xv,fv,xv,fpls,'--',(xv+0.5*dx),RGHTFLX,'o','LineWidth',2)
%   grid on
%   legend('$f(u)$','$f^+(u)$','$f^+_{WENO}$','Location','best')
%   xlabel('$x$')
%   ylabel('$f(u)$')
%   db=1;
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% ------- FLUX FOR LEFT RUNNING WAVES ----- %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% WENO using appropriate stencils for the flux from the left running wave, f-

% Linear weights for O(dx^3) stencil
d0=3/10; % ENSURE LINEAR WEIGHTS ARE CORRECT FOR f- TERMS
d1=3/5;
d2=1/10;

fm1=fmns(2);
fc0=fmns(3);
fp1=fmns(4);
fp2=fmns(5);
fp3=fmns(6);


% Weights for ENO Stencils
c0j = [1/3 5/6 -1/6]; % I+1,I,I-1 stencil (r=0 stencil)
c1j = [-1/6 5/6 1/3]; % I+2,I+1,I stencil (r=1 stencil)
c2j = [1/3 -7/6 11/6]; % I+3,I+2,I+1 stencil (r=2 stencil)

S0 = [fp1 fc0 fm1];
S1 = [fp2 fp1 fc0];
S2 = [fp3 fp2 fp1];



fr0 = c0j*S0';
fr1 = c1j*S1';
fr2 = c2j*S2';

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

LEFTFLX = w0*fr0 + w1*fr1 + w2*fr2;


FM1=RGHTFLX+LEFTFLX;


%%% Visualize RECONSTRUCTED!!! split flux if needed
% if iplt==1
%   figure(9)
%   plot(xv,fv,xv,fmns,'--',(xv+0.5*dx),LEFTFLX,'o','LineWidth',2)
%   grid on
%   legend('$f(u)$','$f^+(u)$','$f^-_{WENO}$','Location','best')
%   xlabel('$x$')
%   ylabel('$f(u)$')
%   db=1;
% 
% 
%   figure(10)
%   plot(xv,fv,xv,fmns,'--',xv,fpls,'--',(xv+0.5*dx),FLUXFACE,'o','LineWidth',2)
%   grid on
%   legend('$f(u)$','$f^-(u)$','$f^+$','$\hat{f}_{j+{1/2}}$','Location','best')
%   xlabel('$x$')
%   ylabel('$f(u)$')
%   db=1;
% end

% FP1=FLUXFACE;
% FM1=circshift(FP1,1);
% FM1(1)=FM1(il); % PERIODIC BC HARD FIX

RHSLOC=FP1-FM1;
db=1;

end