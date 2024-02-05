function Earth3d(Re)

% Makes 3D representation of the Earth.
% 
% PROTOTYPE
%   Earth3d(Re)
% 
% INPUT:
%   Re = radius of the Earth [km]
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

% If not otherwise specified, the radius is 6378 km by default:

if nargin == 0
    Re=6378;
end

I = imread('EarthTex.jpg');
I = imrotate(I,180);
I = flip(I,2);

axis equal

hold on
image([0 360],[-90 90],I,'CDataMapping', 'scaled');
[x,y,z] = sphere(100);

props.FaceColor= 'texture';
props.EdgeColor = 'none';
props.Cdata = I;

surface(x*Re,y*Re,z*Re,props);