function [rvect,vvect,p,r,vpf] = par2car(a,e,i,Omega,omega,theta,mu)

% Converts Keplerian parameters into Cartesian coordinates and velocities.
% 
% PROTOTYPE:
%   [rvect,vvect,p,r,vpf] = par2car(a,e,i,Omega,omega,theta,mu)
% 
% INPUT:
%   a[1] = semi-major axis [km]
%   e[1] = eccentricity [-]
%   i[1] = inclination [rad]
%   Omega[1] = RAAN [rad]
%   omega[1] = argument of periapsis [rad]
%   theta[1] = true anomaly [rad]
%   mu[1] = gravitational parameter [km^3/s^2]
%   
% (alternatively, if nargin <3)
%   a[6] = Keplerian parameters set
%   e[1] = gravitational parameter [km^3/s^2]
% 
% OUTPUT:
%   rvect[3] = position vector, can either be row or column [km]
%   vvect[3] = velocity vector, can either be row or column [km/s]
%   p[1] = semi-latus rectum [km]
%   r[1] = scalar distance from attractor [km]
%   vpf[3] = velocity vector (column) in perifocal reference frame [km/s]
% 
% FUNCTIONS CALLED:
%   (none)
%
% CONTRIBUTORS:
%    Luca Aufiero, Andrea Caushi, Matteo Luciardello Leccardi
% 
% -------------------------------------------------------------------------

% If the user doesn't define mu in the inputs, default to Earth's value:

if nargin < 7 
    mu = 398600;
end

% If the user inputs Keplerian parameters as a vector, as opposed to single
% scalar values:

if nargin < 3
    if nargin == 1
        mu = 398600;
    else
        mu = e;
    end
    % a = a(1);
    e = a(2);
    i = a(3);
    Omega = a(4);
    omega = a(5);
    theta = a(6);
    a = a(1);
end

% if abs(e) > 1
%     warning('abs(e) > 1')
% end

p = a*(1-e^2);
r = p/(1+e*cos(theta));
rpf = r*[cos(theta);sin(theta);0];
vpf = sqrt(mu/(p))*[-sin(theta);e+cos(theta);0];

ROmega = [cos(Omega) sin(Omega) 0
         -sin(Omega) cos(Omega) 0
         0 0 1];

Romega = [cos(omega) sin(omega) 0
         -sin(omega) cos(omega) 0
         0 0 1];

Ri = [1 0 0
      0 cos(i) sin(i)
      0 -sin(i) cos(i)];

T = Romega*Ri*ROmega;

rvect = T'*rpf;
vvect = T'*vpf;
