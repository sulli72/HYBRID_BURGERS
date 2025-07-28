function UPLT(u,xv,ishk)

sqpos=[50 50 800 700];

f=figure(1);
plot(xv,u,'k-o','LineWidth',1.5) %could be either point values or cell avgs 
grid on

% ishk=find(psi==1);

hold on
% plot(xv(ishk),u(ishk),'r-','LineWidth',1.5)
plot(xv(ishk),u(ishk),'ro','MarkerFaceColor','r','MarkerEdgeColor','r','MarkerSize',5)
hold off

xlabel('$x$')
ylabel('$u(x)$')
xlim([-pi pi]);
ylim([ -0.6 1.6]);
% ylim([-3.5 3.5]);
xticks([-pi -pi/2 0 pi/2 pi]);
xticklabels({'$-\pi$' '$-\pi/2$' '$0$' '$\pi/2$' '$\pi$'})
FORMATFIG(f,sqpos,0,0,0);


end