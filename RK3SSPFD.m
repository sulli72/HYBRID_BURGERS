function UNP1 = RK3SSPFD(un,xv,dt,psi,xix,n)
% 3rd Order Runge-Kutta Explicit Time Marching
% Inputs:
% un == solution at current timestep
% xv == grid vector
% dt == timestep size
% psi == vector of shock detection sensor values
% xix == curvlinear metrics for stretched grids

% Outputs:
% UNP1 == solution after RK3 algorithm marches forward in time

%%% General Curvilinear Conservative formulation: 
%  UNP1(j) = UN(j) - dt/dxi * xi_x ( F(j+1/2) - F(j-1/2) ) 

il=length(xv);
RHS=zeros(il,1);

% Reorient metric vector
xix=xix';


% NOTE: UNCOMMENT FLUX SCHEME THAT YOU WANT TO USE
% (SHOULD USE SAME SCHEME AT EACH RUNGE KUTTA STAGE)

% [RHS] = XWENOFLUX(un,xv);    % <-- 5th ORDER CHARCTERISTIC WENO
% [RHS] = COMPACTDIFF6(un,xv); % <-- 6th ORDER COMPACT DIFFERENFCE
% [RHS] = CENTRAL6(un,xv);     % <-- 6th ORDER CENTRAL
% [RHS] = CENTRAL4(un,xv);     % <-- 4th ORDER CENTRAL
% [RHS] = CENTRAL2(un,xv);     % <-- 2nd ORDER CENTRAL
% [RHS] = UPWIND1(un,xv);      % <-- 1st ORDER UPWIND
% RHS = HYBRIDFLUX(un,psi,xv); % <-- INEFFICIENT HYBRID VERSION
RHS = HYBRIDMASK(un,psi,xv,n);   % <-- MORE EFFICIENT HYBRID VERSION USING POINT MASKING

% Metric Step
RHS = xix.*RHS;

u1 = un - dt*(RHS);

% [RHS] = XWENOFLUX(u1,xv);
% [RHS] = COMPACTDIFF6(u1,xv);
% [RHS] = CENTRAL6(u1,xv);
% [RHS] = CENTRAL4(u1,xv);
% [RHS] = CENTRAL2(u1,xv);
% [RHS] = UPWIND1(u1,xv);
% RHS = HYBRIDFLUX(u1,psi,xv);
RHS = HYBRIDMASK(u1,psi,xv,n);

% Metric Step
RHS = xix.*RHS;

u2 = 3/4*un + 1/4*(u1 - dt*(RHS));

% [RHS] = XWENOFLUX(u2,xv);
% [RHS] = COMPACTDIFF6(u2,xv);
% [RHS] = CENTRAL6(u2,xv);
% [RHS] = CENTRAL4(u2,xv);
% [RHS] = CENTRAL2(u2,xv);
% [RHS] = UPWIND1(u2,xv);
% RHS = HYBRIDFLUX(u2,psi,xv);
RHS = HYBRIDMASK(u2,psi,xv,n);


% Metric Step
RHS = xix.*RHS;

UNP1 = 1/3*(un) + 2/3*(u2 - dt*(RHS));


end