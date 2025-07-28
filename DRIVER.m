%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                        %
%                                        %
%  HYBRIDIZED FINITE DIFFERENCE CODE FOR %
%   BURGERS EQUATION ON STRETCHED GRIDS  %
%                                        %
%    ------ By Jack Sullivan ------      %
%                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% SCHEMES:

% Spatial schemes
% 2nd order central
% 4th order central
% 6th order central
% 6th order compact
% 5th order WENO
% Hybrid Central 6-WENO5 
%
% Time scheme
% SSP-RK3 



%% Variable Clearing and Path Appendments
clearvars; close all; clc;
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex');
addpath 'C:\Users\jacks\Documents\MATLAB\Add-Ons\Collections\linspecer';
addpath 'C:\Users\jacks\Documents\MATLAB\Add-Ons\Collections\cmap\cmap-master\cmap-master\';
addpath 'C:\Users\jacks\OneDrive - The Ohio State University\MATLAB_FORMATTING\';
addpath 'C:\Users\jacks\OneDrive - The Ohio State University\SWPISO\3DFLOW_VIS';
sqpos=[50 50 800 700];
nfg=0;


%% Grid Generation
ilvect=[10 20 40 80 160 320]; % <-- grids of various refinement
ngrds=length(ilvect);

% Build legend for plotting later
for igrid=1:length(ilvect)
  il=ilvect(igrid);
  legbase='N =';
  legext=num2str(il);
  legadd=strcat(legbase,legext);
  legvec{igrid}=legadd;
end

%%% ----------- Case parameters -------------- %%%%
% Uncomment desired line to run cases 

% Shock capturing on fine grid with movie
imvie=1; ierror=0; norms=0; tf=2; grd1=5; grd2=grd1; dtflag=1; svflg=0;

% Convergence of smooth solutions over all desired grids
% imvie=0; ierror=1; norms=1; tf=.5; grd1=1; grd2=ngrds; dtflag=0; svflg=0;


%%% ----------- CODE BEGIN ------------------ %%%%

% Grid loop -- considered for convergence studies
for igrid=grd1:grd2

  %%%% ------ PROBLEM SETUP ------- %%%%
  x0=-pi; xf=pi;

  %%% Create IC
  il=ilvect(igrid);
  a0=0.5;  a1=1;

  %%% GRID GENERATION
  grdopt = 1; % <-- constant spacing (use for convergence tests)
  % grdopt = 2; % <-- fine spacing near shock, coarse elsewhere
  % grdopt = 3; % <-- coarse spacing near shock, fine elsewhere
  % grdopt = 4; % <-- 'randomized' grid to test metric stability

  [xv] = GRIDMKR(xf,x0,il,grdopt);

  %%% METRIC COMPUTATION
  [xix] = METRIC(xv);

  %%% INITIAL CONDITION MAKER
  [u0]=FDICMKR(xv,a0,a1);


  % COMPUTE LINEAR AND QUADRATIC INVARIANTS OF INITIAL CONDITION
  UOFT(1) = trapz(xv,u0);
  EOFT(1) = trapz(xv,0.5*(u0.^2));


  %%% Time Marching Parameters
  CFL=.5;
  dxv = xv(2:end)-xv(1:end-1);
  dt=CFL*min(dxv);
  if dtflag==0
    dx=min(dxv);
    dt=max(abs(u0))*dx*dx;
  end
  tv=0:dt:tf;
  nsteps=length(tv);

  % Fix time vector for exact solution comparison
  if tv(end)<tf
    dtfin=tf-tv(end);
    tv(end+1)=tv(end)+dtfin;
    nsteps=nsteps+1;
  end

  % Apply IC
  un=u0; un=un';
  USOLN=zeros(il,nsteps-1);
  ISHK=zeros(il,nsteps-1);
  %%%% ------ MAIN CODE LOOP ------- %%%%
  for n=1:nsteps-1

    % Compute local dt
    mydt=tv(n+1)-tv(n);

    % Call shock sensor
    [ishk]=SENSOR(un); % <-- EDIT SHOCK CAPTURING PARAMETERS IN HERE

    % Add new solution to field
    USOLN(:,n)=un;
    ISHK(ishk,n)=1;

    % Call Time Marching scheme
    unp1=RK3SSPFD(un,xv,mydt,ishk,xix,n); % <-- SELECT SPATIAL SCHEMES IN HERE
    un=unp1; % <-- update solution


    % COMPUTE LINEAR AND QUADRATIC INVARIANTS - FOR CONSERVATION CHECKS
    UOFT(n+1) = trapz(xv,un);
    EOFT(n+1) = trapz(xv,0.5*(un.^2));

    % Plotting if desired
    if imvie==1
      UPLT(un,xv,ishk);
      pause(0.005)
    end

  end
  %%%% ------ END MAIN CODE LOOP ------- %%%%

  %%%% ------ ALL BELOW IS POST-PROCESSING ------- %%%%

  %% Check conservation of linear and quadratic invariants 

  f=figure(2);
  plot(tv,UOFT/UOFT(1),'r-',tv,EOFT/EOFT(1),'b--','LineWidth',2)
  xlabel('$t$')
  ylabel('$I(u)$')
  legend('$\int_\Omega u(t)$','$\int_\Omega |u(t)|_2^2$','Location','best');
  title('Evolution of Burgers Equation Invariants')
  ylim([0 2])
  FORMATFIG(f,sqpos,0,0,0)


  %% Convergence Calculation - Smooth Solutions Only
  if ierror==1
  
    % Get approximate soln
    unum=un;
    tend=tv(nsteps);

    % Determine exact solution, error in various norms
    [linf,l1,l2,uex]=FDERRCOMP(xv,tend,unum,a0,a1);


    % Compare numerical and exact solutions at each grid point
    iexplt=1;
    if iexplt==1
      f=figure(3);
      plot(xv,unum,'o',xv,uex,'s','LineWidth',1)
      xlim([-pi pi]);
      xticks([-pi -pi/2 0 pi/2 pi]);
      xticklabels({'$-\pi$' '$-\pi/2$' '$0$' '$\pi/2$' '$\pi$'})
      legend('Numerical','Exact','Location','best');
      FORMATFIG(f,sqpos,0,0,0)

      ferr=figure(4);
      einf=abs(unum-uex);
      semilogy(xv,einf,'-','LineWidth',1.2)
      title('Absolute Error')
      hold on
      xlim([-pi pi]);
      xticks([-pi -pi/2 0 pi/2 pi]);
      xticklabels({'$-\pi$' '$-\pi/2$' '$0$' '$\pi/2$' '$\pi$'})
      ylim([1e-16 2]);
      legend(legvec{1:igrid},'Location','best');

      if igrid==grd2
        FORMATFIG(ferr,sqpos,0,0,0)
      end
    end

    % Get evolution of error over various grids
    if norms==1 && grd1~=grd2
      if igrid==1
        oldl2=l2;
        oldl1=l1;
        oldlinf=linf;
        k1=0;k2=0;kinf=0;
        % fprintf(1,'%i %1.14f %1.14f %f\n',il,linf,l1,k1);

      else
        newl2=l2;
        newl1=l1;
        newlinf=linf;

        erat2=oldl2/newl2;
        erat1=oldl1/newl1;
        eratinf=oldlinf/newlinf;

        k2=log2(erat2);
        k1=log2(erat1);
        kinf=log2(eratinf);
        % fprintf(1,'%i %1.14f %1.14f %f\n',il,linf,l1,k1);
        fprintf(1,'L1 convergence rate = %f\n',k1)

        l1vec(igrid)=newl1;

        oldl2=newl2;
        oldl1=newl1;
        oldinf=newlinf;

      end
    end
  end




end % grid loop

if exist('ferr','var') == 1
  FORMATFIG(ferr,sqpos,0,0,0)
end

%% Order Rate Convergence

% Plot convergence rates of smooth solutions on successively refined grids
if ierror==1
  scl=.12;
  xplot=ilvect(grd1:grd2);
  xref=ilvect;

  y1ref=xref.^-1;  y1pls=(scl*xref).^-1;   y1ref=y1ref+y1pls;
  y2ref=xref.^-2;  y2pls=(scl*xref).^-2;   y2ref=y2ref+y2pls;
  y4ref=xref.^-4;  y4pls=(scl*xref).^-4;   y4ref=y4ref+y4pls;
  y5ref=xref.^-5;  y5pls=(scl*xref).^-5;   y5ref=y5ref+y5pls;
  y6ref=xref.^-6;  y6pls=(scl*xref).^-6;   y6ref=y6ref+y6pls;


  f=figure(5);
  loglog(xref,y1ref,'k-v',xref,y2ref,'k-s',xref,y4ref,'k-o','LineWidth',1)
  hold on
  loglog(xref,y5ref,'k-^',xref,y6ref,'k-d','LineWidth',1)
  loglog(xplot,l1vec,'r-*','LineWidth',1.5)
  hold off
  title('$L1 - Convergence$')
  xlim([ilvect(grd1)-5 ilvect(grd2)+.25*ilvect(grd2)]);
  ylim([1e-10 1e0]);
  % grid on
  xlabel('$N_x$')
  ylabel('$L_1$')
  legend('$O(h)$','$O(h^2)$','$O(h^4)$','$O(h^5)$','$O(h^6)$','Current Scheme','Location','best')
  FORMATFIG(f,sqpos,1,1,0)
end

%% ANIMATE FIELD
imovie=1;

if imovie==1

  nfg=nfg+1;
  f=figure(10);
  % nfnw=ceil(nt/2);
  nt=nsteps-1;
  nfnw=nt;

  for n=1:nfnw
    un=USOLN(:,n);
    ishk=find(ISHK(:,n)==1);


    % UPLT(un,xv,ishk);

    plot(xv,un,'k-o','LineWidth',1.5) %could be either point values or cell avgs
    grid on

    % ishk=find(psi==1);

    hold on
    % plot(xv(ishk),u(ishk),'r-','LineWidth',1.5)
    plot(xv(ishk),un(ishk),'ro','MarkerFaceColor','r','MarkerEdgeColor','r','MarkerSize',5)
    hold off

    xlabel('$x$')
    ylabel('$u(x)$')
    xlim([-pi pi]);
    ylim([ -0.6 1.6]);
    % ylim([-3.5 3.5]);
    xticks([-pi -pi/2 0 pi/2 pi]);
    xticklabels({'$-\pi$' '$-\pi/2$' '$0$' '$\pi/2$' '$\pi$'})
    fontsize(gcf,16,'points')
    FORMATFIG(f,sqpos,0,0,0);


    drawnow
    frame = getframe(f);
    im{n} = frame2im(frame);

    % fontsize(gcf,scale=1.3)

    % exportgraphics(gcf,'test.gif','Append',true);
    fprintf(1,'DID FRAME WRITE %i/%i\n',n,nfnw)
  end
end

%% GIF WRITING
% filename = "blowup_sensor_evolution.gif"; % Specify the output file name
dtime = 0.1; % <-- delay time in seconds
igif=1;

if igif==1
  for idx = 1:nfnw
    [A,map] = rgb2ind(im{idx},256);
    if idx == 1
      imwrite(A,map,filename,"gif",LoopCount=Inf, ...
        DelayTime=dtime)
    else
      imwrite(A,map,filename,"gif",WriteMode="append", ...
        DelayTime=dtime)
    end
    fprintf(1,'DID GIF FRAME %i/%i\n',idx,nfnw)
  end


end