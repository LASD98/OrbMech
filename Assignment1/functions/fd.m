function [Fd,d,e_minus,e_plus,d_minus,d_plus] = fd(rp,vinf_minus,vinf_plus,mu)
% rp solver function for fzero or fsolve, for hyperbolic powered gravity assist
% INPUTs:
%         rp[1]: unknown pericenter radius [km]
%         vinf_minus,vinf_plus[3x1]: incoming & outcoming velocities 
%         mu[1]: planetary constant [km^3/s^2]
% OUTPUTs:
%         Fd: equation to be solved
%         d: total turning angle [deg] 
%         e_minus, e_plus: eccentricity of each arc of hyperbola
%         d_minus, d_plus: turning angle provided by each arc

e_minus = 1+rp*norm(vinf_minus)^2/mu;
d_minus = 2*asind(1/e_minus);

e_plus = 1+rp*norm(vinf_plus)^2/mu;
d_plus = 2*asind(1/e_plus);

d = acosd(dot(vinf_minus,vinf_plus)/(norm(vinf_minus)*norm(vinf_plus)));

Fd = d - (d_minus+d_plus)/2;

end