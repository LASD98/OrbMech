function [p,pp] = plotOrbit(a,e,i,Omega,omega,theta0,thetaf,dtheta)

% Plots orbits, given Keplerian parameters, true anomaly boundaries and
% discretization step.
% 
% PROTOTYPE:
%   [p,pp] = plotOrbit(a,e,i,Omega,omega,theta0,thetaf,dtheta)
% 
% INPUT:
%   a[1] = semi-major axis [km]
%   e[1] = eccentricity [-]
%   i[1] = inclination [rad]
%   Omega[1] = RAAN [rad]
%   omega[1] = argument of periapsis [rad]
%   theta0[1] = lower true anomaly boundary [rad]
%   thetaf[1] = upper true anomaly boundary [rad]
%   dtheta[1] = true anomaly discretization step [rad]
% 
% OUTPUT:
%   p[1] = graphical object associated with the plotted orbit (line element)
%   pp[1] = graphical object associated with the apse line (line element)
% 
% FUNCTIONS CALLED:
%   par2car
% 
% CONTRIBUTORS:
%   Luca Aufiero, Andrea Caushi, Matteo Luciardello Leccardi
% 
% -------------------------------------------------------------------------

theta = theta0:dtheta:thetaf;

n = size(theta,2);
state = zeros(n,7);

for j = 1:n
    
    thetaj = theta(j);
    [r,v] = par2car(a,e,i,Omega,omega,thetaj);
    state(j,:) = [thetaj r' v'];
    
end

hold on
p = plot3(state(:,2),state(:,3),state(:,4),'LineWidth',0.5);

% [rr1,~] = par2car(a,e,i,Omega,omega,0);
% [rr2,~] = par2car(a,e,i,Omega,omega,pi);

% pp = plot3([rr1(1),rr2(1)],[rr1(2),rr2(2)],[rr1(3),rr2(3)],'LineWidth',0.5,'LineStyle','--','Color','k');
% pp = [];
% p = fill3(state(:,2),state(:,3),state(:,4),[0 1 1]);
% set(p,'facealpha',.3)

% axis vis3d
% hold off