function dy = orbitJ2(~,y,mu,Re,J2)

% Defines the two-body differential problem with J2 perturbations, to be
% solved with an ODE solver, like ode113.
% 
% PROTOTYPE
%   dy = orbitJ2(~,y,mu,Re,J2)
% 
% INPUT:
%   y[6] = vector of Cartesian position and velocities, positions first.
%   mu[1] = gravitational parameter [km^3/s^2]
%   Re[1] = radius of the parent body [km]
%   J2[1] = second zonal harmonic of the zonal variations of the gravitational
%     field, due to the oblateness of the attractor [-]
% 
% OUTPUT:
%   dy[6] = vector of derivatives of Cartesian positions and velocities,
%     positions first.
% 
% FUNCTIONS CALLED:
%   (none)
% 
% CONTRIBUTORS:
%   Luca Aufiero, Andrea Caushi, Matteo Luciardello Leccardi
% 
% -------------------------------------------------------------------------

% If the user doesn't define Re and J2 in the inputs, default to Earth's value:

if nargin < 4
    Re = 6378.137;
    J2 = 0.00108263;
end

% If the user doesn't define mu in the inputs, default to Earth's value:

if nargin < 3
    mu = 398600;
end

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

dy = [
      y(4)
      y(5)
      y(6)
      -mu*y(1)/r^3 + aJ2(1)
      -mu*y(2)/r^3 + aJ2(2)
      -mu*y(3)/r^3 + aJ2(3)
      ];