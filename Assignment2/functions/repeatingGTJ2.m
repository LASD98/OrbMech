function aJ2 = repeatingGTJ2(m,k,a,e,i,omE,J2,Re,mu)

% Computes semi-major axis required to get a repeating ground track every k
% satellite revolutions and m parent body revolutions.
% 
% PROTOTYPE:
%   aJ2 = repeatingGTJ2(m,k,a,e,i,omE,J2,Re,mu)
% 
% INPUT:
%   m[1] = required parent body revolutions
%   k[1] = required satellite revolutions
%   a[1] = semi-major axis of generic configuration [km] 
%   e[1] = eccentricity
%   i[1] = inclination [rad]
%   omE[1] = parent body rotation speed [rad/s]
%   J2[1] = second zonal harmonic of the zonal variations of the gravitational
%     field, due to the oblateness of the attractor [-]
%   Re[1] = radius of parent body [km]
%   mu[1] = gravitational parameter [km^3/s^2]
% 
% OUTPUT:
%   aJ2[1] = semi-major axis that makes ground tracks repeat [km]
% 
% FUNCTIONS CALLED:
%   (none)
% 
% CONTRIBUTORS:
%   Luca Aufiero, Andrea Caushi, Matteo Luciardello Leccardi
% 
% -------------------------------------------------------------------------

if nargin < 9
    mu = 398600; % default is Earth's gravitational parameter
    if nargin < 8
        Re = 6378; % default is Earth's radius
        if nargin < 7
            J2 = 0.00108263; % default is Earth's second zonal harmonic
            if nargin < 6
                omE = deg2rad(15.04/3600); %  default is Earth's rotation speed
            end
        end
    end
end


n = @(a) sqrt(mu/a^3);
    
M0dot = @(a) -3/2*(sqrt(mu)*J2*Re^2)/((1-e^2)^2*a^(7/2))*(1-3/2*sin(i)^2);
    
Omdot = @(a) -3/2*(sqrt(mu)*J2*Re^2)/((1-e^2)^2*a^(7/2))*cos(i);
    
omdot = @(a) -3/2*(sqrt(mu)*J2*Re^2)/((1-e^2)^2*a^(7/2))*(5/2*sin(i)^2-2);
    
    
fun = @(a) m/k - (omE - Omdot(a))/(n(a)+omdot(a)+M0dot(a));
    
    
aJ2 = fzero(fun,a); % generic semi-major axis is used as initial guess