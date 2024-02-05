function dy = orbitSRPJ2(t,y,mu,Re,J2,omE,body)

% Defines the two-body differential problem with SRP and J2 perturbations,
% to be solved with an ODE solver, like ode113.
% 
% PROTOTYPE
%   dy = orbitSRPJ2(t,y,mu,Re,J2,omE,body)
% 
% INPUT:
%   y[6] = vector of Cartesian position and velocities, positions first.
%   mu[1] = gravitational parameter [km^3/s^2]
%   Re[1] = radius of the parent body [km]
%   J2[1] = second zonal harmonic of the zonal variations of the gravitational
%     field, due to the oblateness of the attractor [-]
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

aJ2 = [3/2*J2*mu*Re^2/r^4*(xx/r*(5*zz^2/r^2-1))
       3/2*J2*mu*Re^2/r^4*(yy/r*(5*zz^2/r^2-1))
       3/2*J2*mu*Re^2/r^4*(zz/r*(5*zz^2/r^2-3))];

[RR1, ~] = kep2car(kep(1), kep(2), 23.4*pi/180, kep(4), kep(5), kep(6), mu_s);
rrSSC = RR1(:) + [xx;yy;zz];

rSSC = norm(rrSSC);
aSRP = pSR*(AU^2/norm(rSSC)^2)*cr*(A_m)*1e-3;
a_SRP_x = -aSRP*rrSSC(1)/norm(rSSC);
a_SRP_y = -aSRP*rrSSC(2)/norm(rSSC);
a_SRP_z = -aSRP*rrSSC(3)/norm(rSSC);

dy = [
      y(4)
      y(5)
      y(6)
      -mu*y(1)/r^3 + a_SRP_x + aJ2(1)
      -mu*y(2)/r^3 + a_SRP_y + aJ2(2)
      -mu*y(3)/r^3 + a_SRP_z + aJ2(3)
      ];