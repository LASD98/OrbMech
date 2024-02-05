function v_rot = rodrigues(v,u,delta)

%--------------------------------------------------------------------------
%
% Rotates a vector in a certain plane of the angle delta
%
% PROTOTYPE:
% v_rot = rodrigues(v,u,delta)
%
% INPUT:
% v     [1] = Starting vector
% u     [1] = Direction normal to the rotation plane
% delta [1] = Angle of rotation
%
% OUTPUT:
% v_rot [1] = Rotated vector

%--------------------------------------------------------------------------

v_rot = (v.*cos(delta)+cross(u,v).*sin(delta)+u.*(dot(u,v)*(1-cos(delta))));

end