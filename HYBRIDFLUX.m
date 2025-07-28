function [RHS] = HYBRIDFLUX(utmp,psi,xv)

[RHSU] = XWENOFLUX(utmp,xv);
[RHSC] = CENTRAL6(utmp,xv);

RHS = (1-psi).*RHSC + psi.*(RHSU);

end