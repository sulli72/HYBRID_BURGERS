function FLXFNC = U2FLX(u)
  FLXFNC = 0.5*(u.^2); % <-- non-linear flux
  % FLXFNC = 1*u; % <-- linear flux
end