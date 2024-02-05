function [a,e,i,Om,om,theta] = car2kep(rr,vv,mu)

% Converts Cartesian coordinates and velocities into Keplerian parameters.
% 
% PROTOTYPE:
%   [a,e,i,Om,om,theta] = car2kep(rr,vv,mu)
% 
% INPUT:
%   rr[3] = position vector, can either be row or column [km]
%   vv[3] = velocity vector, can either be row or column [km/s]
%   mu[1] = gravitational parameter [km^3/s^2]
% 
% OUTPUT:
%   a[1] = semi-major axis [km]
%   e[1] = eccentricity [-]
%   i[1] = inclination [rad]
%   Om[1] = RAAN [rad]
%   om[1] = argument of periapsis [rad]
%   theta[1] = true anomaly [rad]
% 
% FUNCTIONS CALLED:
%   (none)
% 
% CONTRIBUTORS:
%    Luca Aufiero, Andrea Caushi, Matteo Luciardello Leccardi
% 
% -------------------------------------------------------------------------

r = norm(rr);
v = norm(vv);
k = [0 0 1]';

a = -mu / (v^2 - 2*mu/r);

hh = cross(rr,vv);
h = norm(hh);

ee = cross(vv,hh)/mu - rr/r;
e = norm(ee);

i = acos(hh(3)/h);

N = cross(k,hh) / norm(cross(hh,k));

Om = acos(N(1));
if N(2) < 0
    Om = 2*pi - Om;
end

om = acos(dot(N,ee)/e);
if ee(3) < 0
    om = 2*pi - om;
end

theta = acos(dot(rr,ee)/(r*e));
if dot(rr,vv) < 0
    theta = 2*pi - theta;
end




end

