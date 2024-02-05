function [f,Y,Fs] = spectralAnalysis(var,time,LineWidth,flag)

% Performs the spectral analysis of the input signal var, yielding amplitude.
% 
% PROTOTYPE
%   [f,Y,Fs] = spectralAnalysis(var,time,flag)
% 
% INPUT:
%   var[n] = generic signal to analyze [-]
%   time[n] = time of measurement vector [s]
%   flag[1] = defines whether the plot of the signal in the frequency
%       domain is shown or not. 
% 
% OUTPUT:
%   f[n] = frequency range of domain [s^-1]
%   Y[n] = amplitude of signal as a function of frequency [-]
%   Fs[1] = sampling frequency [s^-1]
% 
% FUNCTIONS CALLED:
%   (none)
% 
% CONTRIBUTORS:
%   Luca Aufiero, Andrea Caushi, Matteo Luciardello Leccardi
% 
% -------------------------------------------------------------------------

if nargin < 4
    flag = 1;
    if nargin < 3
        LineWidth = 0.5;
    end
end

mu = 398600;

N = length(var); % tweak this value here

dt = time(end)-time(1); % tweak this value here

T = dt/(N-1);

Fs = 1/T;

df = 1/dt;

yy = var;

y = fft(yy);

ynorm = abs(y/N);

Y = fftshift(ynorm);

f = (-N/2:N/2-1)*(Fs/N);

if flag
    
    semilogy(f,Y,'LineWidth',LineWidth)
    xlabel('Frequency [Hz]')
    ylabel('Amplitude')

    grid on
end
%%
% % This shows that peaks correspond to multiples of the period.
% 
% period = 2*pi*sqrt(var(1)^3/mu);
% 
% % hold on
% for ii = 1:23
%     
%     plot([ii/period,ii/period],[1e-4,1e4],'Color','r')
%     
% end