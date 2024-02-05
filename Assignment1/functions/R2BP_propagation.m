function [t,x] = R2BP_propagation(x0,tspan,mu,J2)
% Restricted 2 body problem propagation
% INPUT:
%   x0[6x1] Initial State of the body ( rx, ry, rz, vx, vy, vz) [ L, L/T ]
%   tspan[N] Propagation time, boundaries[2] or vector[N] can be borth provided [s]
%   mu[1] Gravitational constant for celestial bodies [km^3/s^2]
%   J2[1]  boolean value, adds J2 perturbation
%
% OUTPUT
%   t[N] Propagation time vector [ s ]
%   x[N,6] Propagated state
% 
% FUNCTIONS CALLED:
% ode_2bp.m

if nargin < 4
    J2 = 0;
end
if nargin < 3
    mu = astroConstants(13); % Earth's gravitational parameter [km^3/s^2]
end

options = odeset('RelTol',1e-13,'AbsTol',1e-14);
[t,x] = ode113(@ode_2bp,tspan,x0,options,J2,mu);
end