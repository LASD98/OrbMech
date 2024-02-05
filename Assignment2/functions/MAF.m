function y = MAF(x,N)

% Filters the input signal with a symmetric Moving Average Filter, picking
% N/2-1 points to the left and N/2-1 points to the right of any given point
% to compute the average.
% 
% PROTOTYPE
%   y = MAF(x,N)
% 
% INPUT:
%   x[n] = signal to filter
%   N[1] = number of points required for the average
% 
% OUTPUT:
%   y[n] = filtered signal
% 
% FUNCTIONS CALLED:
%   filtfilt
% 
% CONTRIBUTORS:
%   Luca Aufiero, Andrea Caushi, Matteo Luciardello Leccardi
% 
% -------------------------------------------------------------------------

if nargin < 2
    N = 5000;
end

v = ones(1,round(N/2));
B = conv(v, v);
A = sum(B);

y = filtfilt(B,A,x);
