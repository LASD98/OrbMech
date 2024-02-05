function globe(name,R,r_input)

% plots globe (sphere) with specified texture
% INPUTS:
%       name[char/string]: 'name.extension' of texture file
%       R[1]: radius of globe [km]
%       r_input[3]: position of center of globe [km], optional

if nargin < 3
    r_input = zeros(3,1);
end
% R = astroConstants(23);
[X,Y,Z] = sphere(100);
I = imread(name); % I = flip(I,1);
props.FaceColor = 'texture';
props.EdgeColor = 'none';
props.CData = I;
surf(-R*X+r_input(1),-R*Y+r_input(2),-R*Z+r_input(3),props), hold on, grid on
axis equal padded, pbaspect([1,1,1])

end
