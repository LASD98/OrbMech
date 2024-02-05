function planisphere()

% Represents a map of the surface of Earth.
% 
% PROTOTYPE:
%   planisphere()
% 
% INPUT:
%   (none)
% 
% OUTPUT:
%   (none)
% 
% FUNCTIONS CALLED:
%   (none)
% 
% CONTRIBUTORS:
%   Luca Aufiero, Andrea Caushi, Matteo Luciardello Leccardi
% 
% -------------------------------------------------------------------------

figure
I = imread("EarthTex.jpg");
I = flip(I,1);

image([-180 180],[-90 90],I);
grid on
xticks(-180:30:180)
yticks(-90:30:90)
set(gca,'YDir','normal')
hold on

xlabel('Longitude [deg]')
ylabel('Latitude [deg]')