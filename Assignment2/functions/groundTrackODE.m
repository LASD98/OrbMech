function [alpha, delta, lon, lat, Y] = groundTrackODE(state, lonG0, time, omE, mu)

% Without considering J2 effects, computes right ascension, declination,
% latitude, and longitude, given Keplerian parameters and elapsed time.
% 
% PROTOTYPE:
%   [alpha, delta, lon, lat, Y] = groundTrackODE(state, lonG0, time, omE, mu)
% 
% INPUT:
%   state[6] = Keplerian parameters state vector
%   lonG0[1] = longitude of Greenwich meridian at initial time in inertial
%       reference [rad]
%   time[n] = time vector [s]
%   omE[1] = parent body rotation speed [rad/s]
%   mu[1] = gravitational parameter [km^3/s^2]
% 
% OUTPUT:
%   alpha[n] = right ascension vector [rad]
%   delta[n] = declination vector [rad]
%   lon[n] = longitude vector [rad]
%   lat[n] = latitude vector [rad]
%   Y[nx6] = Cartesian coordinates vector
% 
% FUNCTIONS CALLED:
%   par2car
%   orbit
% 
% CONTRIBUTORS:
%   Luca Aufiero, Andrea Caushi, Matteo Luciardello Leccardi
% 
% -------------------------------------------------------------------------

if nargin < 5 
    mu = 398600; % default is Earth's gravitational parameter
	if nargin < 4
        omE = deg2rad(15.04/3600); % default is Earth's rotation speed
	end
end

a = state(1);
e = state(2);
i = state(3);
Om = state(4);
om = state(5);
th = state(6);

[rr,vv] = par2car(a,e,i,Om,om,th);

opts = odeset('Reltol',1e-12,'AbsTol',1e-15);

[T,Y] = ode113(@(t,y) orbit(t,y,mu), time, [rr;vv], opts);

alpha = zeros(length(time),1);
delta = zeros(length(time),1);
lon = zeros(length(time),1);
lat = zeros(length(time),1);
lonG = zeros(length(time),1);

for jj = 1:length(time)
    
    lonG(jj) = mod(lonG0 + (T(jj)-T(1))*omE, 2*pi);
    
    x = Y(jj,1);
    y = Y(jj,2);
    z = Y(jj,3);
    
    r = norm([x y z]);
    
    delta(jj) = asin(z/r);

    alpha(jj) = acos((x/r)/cos(delta(jj)));

    if y/r < 0
        alpha(jj) = 2*pi-alpha(jj);
    end

    lon(jj) = alpha(jj) - lonG(jj);
    lat(jj) = delta(jj);

    if lon(jj) > pi
        lon(jj)  = lon(jj) - 2*pi;
    end
    
    if lon(jj) < -pi
        lon(jj)  = lon(jj) + 2*pi;
    end

%     lon(jj) = mod(lon(jj),2*pi)-pi;

%     fprintf("%d\n",jj);
end