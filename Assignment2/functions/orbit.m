function dy = orbit(~,y,mu)

% Defines the two-body unperturbed differential problem, to be solved with
% an ODE solver, like ode113.
% 
% PROTOTYPE:
%   dy = orbit(~,y,mu)
% 
% INPUT:
%   y[6] = vector of Cartesian position and velocities, positions first.
%   mu[1] = gravitational parameter [km^3/s^2]
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

dy = [
      y(4)
      y(5)
      y(6)
      -mu*y(1)/r^3
      -mu*y(2)/r^3
      -mu*y(3)/r^3
      ];