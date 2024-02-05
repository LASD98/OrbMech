function [a,e,i,OM,om,th] = car2kep (r_vect,v_vect,mu)
%
% Converter: Cartesian Coordinates & Velocities ---> Keplerian Parameters 
%
% DESCRIPTION:
% This code provides a conversion from cartesian coordinates (r_vect) and 
% velocities (v_vect) vectors to Keplerian Parameters of the orbit.
%
%--------------------------------------------------------------------------
% INPUTS:
%   r_vect     [3x1]       Position Vector           [km]
%   v_vect     [3x1]       Velocity Vector           [km/s]
%   mu         [1x1]       Planetary Constant        [km^3][sec-2]
%--------------------------------------------------------------------------
% OUTPUTS:
%   a          [1x1]       Semi-Major Axis           [km]
%   e          [1x1]       Eccentricity              [-]
%   i          [1x1]       Inclination               [rad]
%   OM         [1x1]       RAAN                      [rad]
%   om         [1x1]       Argument of Periapsis     [rad]
%   th         [1x1]       True Anomaly              [rad]
%--------------------------------------------------------------------------
%% Conversion Routine
%Direction and Velocity Norms
r = norm(r_vect);
v = norm(v_vect);
    
%Angular Momentum
h_vect = cross(r_vect, v_vect);
h = norm(h_vect);
%Eccentricity Vector
c1 = cross(v_vect, h_vect);
e_vect = c1./mu - r_vect./r;
e = norm(e_vect);
%Mechanical Energy
E = v^2/2 - mu/r;
%Semi-Major Axis
a = -mu/(2*E);
%Orbit Inclination
k_vect = [0; 0; 1];
c2 = dot(h_vect,k_vect);
i = acos (c2./h);
%Nodes Line Vector
c3 = cross(k_vect,h_vect);
n_vect = c3./norm(c3);
%Right Ascension of the Ascending Node
i_vect = [1;0;0];
j_vect = [0;1;0];
c4 = dot(n_vect,i_vect);        %Projection of Vector n on i
c5 = dot(n_vect,j_vect);        %Projection of Vector n on j
if(c5 >= 0)
    OM = acos(c4/norm(n_vect));
else
    OM = -acos(c4/norm(n_vect)) + 2*pi;
end
%Anomaly of the Pericenter
c6 = dot(e_vect, k_vect);       %Projection of Vector e on k
c7 = dot(n_vect, e_vect);
if (c6>=0)
    om = acos(c7./(e*norm(n_vect)));
else
    om = -acos(c7./(e*norm(n_vect)))+2*pi; 
end
%True Anomaly
c8 = dot(v_vect,r_vect);
c9 = dot(r_vect,e_vect);
if(c8>=0)
    th = acos(c9./(r*e));
else
    th = -acos(c9./(r*e))+2*pi;
end
if nargout == 1
    a = [a,e,i,OM,om,th];
end
end