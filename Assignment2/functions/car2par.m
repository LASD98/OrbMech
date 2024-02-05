function [a,e,i,Om,om,th,hvect,evect] = car2par(rr,vv,mu)

% Converts Cartesian coordinates and velocities into Keplerian parameters.
% 
% PROTOTYPE:
%   [a,e,i,Om,om,th,hvect,evect] = car2par(rr,vv,mu)
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
%   th[1] = true anomaly [rad]
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

% If the user doesn't define mu in the inputs, default to Earth's value:

if nargin < 3
    mu = 398600;
end

a = (2/r - (v^2)/mu)^-1;

hvect = cross(rr,vv);
h = norm(hvect);

evect = cross(vv,hvect)/mu - rr/r;
e = norm(evect);

i = acos(hvect(3)/h);

k = [0,0,1];

if i == 0 % If the orbit is equatorial, to avoid singularities in RAAN:
        
    Om = 0;    
    
    if evect(2) > 0
        om = acos(evect(1)/e);
    else
        om = 2*pi - acos(evect(1)/e);
    end
    
    vr = dot(vv,rr)/r;

    if vr>0
        th = acos(dot(rr,evect)/(r*e));
    else
        th = 2*pi - acos(dot(rr,evect)/(r*e));
    end

elseif e == 0 % If the orbit is circular, to avoid singularities in om:

    Nvect = cross(k,hvect)/norm(cross(k,hvect));

    if Nvect(2)>=0
        Om = acos(Nvect(1));
    else
        Om = 2*pi - acos(Nvect(1));
    end

    om = 0;

    vr = dot(vv,rr)/r;

    if vr>0
        th = acos(dot(rr,Nvect)/r);
    else
        th = 2*pi - acos(dot(rr,Nvect)/r);
    end
    
elseif i == 0 && e == 0 % If the orbit is circular and equatorial:

    Om = 0;
    om = 0;
    
    vr = dot(vv,rr)/r;
    
    if vr>0
        th = acos(rr(1)/r);
    else
        th = 2*pi - acos(rr(1)/r);
    end
    
else % If the orbit is neither circular nor equatorial:

    Nvect = cross(k,hvect)/norm(cross(k,hvect));

    if Nvect(2)>=0
        Om = acos(Nvect(1));
    else
        Om = 2*pi - acos(Nvect(1));
    end

    if evect(3)>=0
        om = acos(dot(Nvect,evect)/e);
    else
        om = 2*pi - acos(dot(Nvect,evect)/e);
    end

    vr = dot(vv,rr)/r;

    if vr>0
        th = acos(dot(rr,evect)/(r*e));
    else
        th = 2*pi - acos(dot(rr,evect)/(r*e));
    end

end