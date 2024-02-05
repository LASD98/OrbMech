function th = anglePostProcessingMod(theta)

% Manipulates a set of angles so that every term belongs to the [0,2pi]
% interval, dividing by 2pi and taking the remainder.
% 
% PROTOTYPE:
%   th = anglePostProcessingMod(theta)
% 
% INPUT:
%   theta[n] = vector of angles to process [-]
% 
% OUTPUT:
%   th[n] = vector of processed angles [-]
% 
% FUNCTIONS CALLED:
%   (none)
% 
% CONTRIBUTORS:
%   Luca Aufiero, Andrea Caushi, Matteo Luciardello Leccardi
% 
% -------------------------------------------------------------------------

for jj = 1:length(theta)
    if theta(jj) > pi
        theta(jj) = mod(theta(jj),2*pi);
    end
%     disp(jj)
end

th = theta;