function om = anglePostProcessingShift(omega)

% Manipulates a set of angles so that every term belongs to the [-2pi,2pi]
% interval, shifting by 2pi.
% 
% PROTOTYPE:
%   om = anglePostProcessingShift(omega)
% 
% INPUT:
%   omega[n] = vector of angles to process [-]
% 
% OUTPUT:
%   om[n] = vector of processed angles [-]
% 
% FUNCTIONS CALLED:
%   (none)
% 
% CONTRIBUTORS:
%   Luca Aufiero, Andrea Caushi, Matteo Luciardello Leccardi
% 
% -------------------------------------------------------------------------

for jj = 1:length(omega)
    if omega(jj) > pi
        omega(jj) = omega(jj) - 2*pi;
    end
%     disp(jj)
end

om = omega;