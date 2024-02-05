function y = ode_2bp(~,x,J2,mu)
%ode_2bp ODE system for the two-body problem (Keplerian motion)
%
% PROTOTYPE
% dy = ode_2bp( t, y )
%
% INPUT:
%   t[1]   Time (can be omitted, as the system is autonomous) [T]
%   y[6x1] State of the body ( rx, ry, rz, vx, vy, vz) [ L, L/T ]
%   J2[1]  boolean value, adds J2 perturbation
%   mu[1] Gravitational parameter of the primary [L^3/T^2]

% OUTPUT:
%   dy[6x1] Derivative of the state [ L/T^2, L/T^3 ]
%
% CONTRIBUTORS:
% 
% VERSIONS:
%

r = x(1:3);
v = x(4:6);
y = zeros(6,1);
y(1:3) = v;

rn = norm(r);
switch J2
    case 0
        aj2 = 0;
    case 1
        Re = astroConstants(23); % mean Earth radius
%       j2 = 0.00108263; % from https://nssdc.gsfc.nasa.gov/planetary/factsheet/earthfact.html
        j2 = astroConstants(9);
        aj2 = 3/2*j2*mu*Re^2/(rn^4)*r/rn.*(5*(r(3)/rn)^2-[1;1;3]); % J2 perturbation
end
y(4:6) = -mu*r/rn^3+aj2; % dynamics
end