function dy = orbitSRP(t,y,mu,Re,omE,body)

% Defines the two-body differential problem with SRP perturbations, to be
% solved with an ODE solver, like ode113.
% 
% PROTOTYPE
%   dy = orbitSRP(t,y,mu,Re,omE,body)
% 
% INPUT:
%   y[6] = vector of Cartesian position and velocities, positions first.
%   mu[1] = gravitational parameter [km^3/s^2]
%   Re[1] = radius of the parent body [km]
%   omE[1] = parent body rotation speed [rad/s]
%   body[2] = vector containing area to mass ratio [m^2/kg] and reflectivity coefficient [-]
% 
% OUTPUT:
%   dy[6] = vector of derivatives of Cartesian positions and velocities,
%     positions first.
% 
% FUNCTIONS CALLED:
%   astroConstants
%   uplanet
%   kep2car
% 
% CONTRIBUTORS:
%   Luca Aufiero, Andrea Caushi, Matteo Luciardello Leccardi
% 
% -------------------------------------------------------------------------

% SRP Constants
AU = astroConstants(2);
pSR = 4.5*1e-6; % SRP @ 1AU
mu_s = astroConstants(4);
A_m = body(1);
cr = body(2);
pl_ID = 3; % Earth
kep = uplanet(t/86400, pl_ID);

xx = y(1);
yy = y(2);
zz = y(3);
% vxx = y(4);
% vy = y(5);
% vzz = y(6);

r = sqrt(xx^2+yy^2+zz^2);

[RR1, ~] = kep2car(kep(1), kep(2), 23.4*pi/180, kep(4), kep(5), kep(6), mu_s);
rrSSC = RR1(:);

rSSC = norm(rrSSC);
aSRP = pSR*(AU^2/rSSC^3)*cr*(A_m)*1e-3;
a_SRP_x = -aSRP*rrSSC(1);
a_SRP_y = -aSRP*rrSSC(2);
a_SRP_z = -aSRP+rrSSC(3);

dy = [
      y(4)
      y(5)
      y(6)
      -mu*y(1)/r^3 + a_SRP_x
      -mu*y(2)/r^3 + a_SRP_y
      -mu*y(3)/r^3 + a_SRP_z
      ];